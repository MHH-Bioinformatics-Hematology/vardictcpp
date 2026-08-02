#include "reference.hpp"
#include "simd.hpp"
#include <stdexcept>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace vardict {

Reference::Reference(const std::string& fastaPath) {
    fai_ = fai_load(fastaPath.c_str());
    if (!fai_) throw std::runtime_error("cannot load faidx for " + fastaPath +
                                        " (is it indexed with `samtools faidx`?)");
}

Reference::~Reference() {
    if (fai_) fai_destroy(fai_);
}

std::string Reference::fetchSeq(const std::string& chr, int& s, int e) const {
    if (s < 1) s = 1;
    int len = 0;
    // faidx_fetch_seq is 0-based, end-inclusive; it clamps e to the contig length.
    char* raw = faidx_fetch_seq(fai_, chr.c_str(), s - 1, e - 1, &len);
    if (!raw || len < 0) { if (raw) free(raw); return std::string(); }
    std::string out(raw, len);
    // Upper-case the loaded window (up to a full chromosome) with SIMD; identical result to the
    // per-byte std::toupper, just 16 bases at a time. Portable across SSE2 / NEON / scalar.
    if (!out.empty()) simd::toupper_ascii(&out[0], out.size());
    free(raw);
    return out;
}

void Reference::fetchWindow(const std::string& chr, int s, int e) {
    if (s < 1) s = 1;
    seq_ = fetchSeq(chr, s, e);   // s is clamped in place
    loadedStart_ = s;
    loadedChr_ = chr;
    contigLen_ = fai_ ? faidx_seq_len(fai_, chr.c_str()) : 0;
    if (contigLen_ < 0) contigLen_ = 0;
    // Java's ReferenceResource.getReference stores referenceSequences only up to sequenceEnd - SEED_1
    // unless the loaded window reaches the contig end (len == sequenceEnd). Reproduce that trailing
    // truncation so has() answers null at the padded window's tail, exactly like ref.get(p) == null.
    int end = loadedEnd();
    primaryEffEnd_ = (contigLen_ > 0 && end >= contigLen_) ? end : end - SEED_1;
    seedBuilt_ = false;
}

void Reference::addSegment(int s, std::string&& seq) {
    if (seq.empty()) return;
    int e = s + (int)seq.size() - 1;
    // Merge with any existing extra segment that overlaps or abuts [s,e] so the list stays a set of
    // disjoint, position-sorted windows (identical bases where they overlap -- both are real reference).
    for (size_t i = 0; i < extra_.size();) {
        Segment& sg = extra_[i];
        int ss = sg.start, se = sg.start + (int)sg.seq.size() - 1;
        if (se + 1 >= s && ss <= e + 1) {           // overlap or adjacency -> fold sg into [s,e]/seq
            int ns = std::min(s, ss), ne = std::max(e, se);
            std::string merged((size_t)(ne - ns + 1), 'N');
            for (int p = ss; p <= se; ++p) merged[p - ns] = sg.seq[p - ss];
            for (int p = s;  p <= e;  ++p) merged[p - ns] = seq[p - s];
            s = ns; e = ne; seq = std::move(merged);
            extra_.erase(extra_.begin() + i);
            continue;
        }
        ++i;
    }
    // Insert keeping extra_ sorted by start.
    size_t pos = 0; while (pos < extra_.size() && extra_[pos].start < s) ++pos;
    extra_.insert(extra_.begin() + pos, Segment{s, std::move(seq)});
}

bool Reference::addGenuine(int a, int b) {
    if (a < 1) a = 1;
    if (b < a) return false;
    if (inGenuine(a, b)) return false;
    genuine_.emplace_back(a, b);
    return true;
}

void Reference::load(const std::string& chr, int start, int end, int pad) {
    int s = start - pad; if (s < 1) s = 1;
    int e = end + pad;
    fetchWindow(chr, s, e);
    // A fresh region load starts a new set of genuine windows and drops any prior extra segments.
    extra_.clear();
    genuine_.clear();
    genuine_.emplace_back(loadedStart_, loadedEnd());
}

void Reference::ensure(int start, int end, int pad) {
    if (loadedChr_.empty()) return;
    // Already covered by the primary window or a previously-loaded extra segment? Then just record the
    // genuinely-requested window (Java would find it via isLoaded and not re-fetch).
    bool covered = (start >= loadedStart_ && end <= loadedEnd()) || inExtra(start, end);
    bool added = addGenuine(start, end);
    if (covered) { if (added) seedBuilt_ = false; return; }
    // Load ONLY a small window around [start,end] as a disjoint segment (Java's getReference window:
    // [start - pad, end + pad], pad = numberNucleotideToExtend + extension), NOT a contiguous gap-fill
    // to the far breakpoint. This is the memory fix: one region reaching a far SV/large-indel breakpoint
    // used to pull in every base in between (multi-Mbp) per worker thread.
    int s = start - pad; if (s < 1) s = 1;
    int e = end + pad;
    addSegment(s, fetchSeq(loadedChr_, s, e));   // s clamped in place; ignored after (segment stores it)
    seedBuilt_ = false;
}

void Reference::buildSeed() const {
    seed17_.clear(); seed12_.clear();
    seedBuilt_ = true;
    if (!extra_.empty()) { buildSeedDisjoint(); return; }
    const int n = (int)seq_.size();
    if (n <= 0) return;
    seed17_.reserve((size_t)n);
    const char* s = seq_.data();

    // Classify all bases to their 3-bit codes once (SIMD), then slide a rolling k-mer encode across the
    // code array. This replaces the previous per-position encodeKmer (which re-read k bases for every
    // position, O(n*k)) with an O(n) pass whose k-mer values are bit-for-bit identical, so the seed
    // maps -- and therefore all downstream calls -- are unchanged. Only positions that could occur more
    // than once are marked 0 (VarDict treats those as non-unique), exactly as before.
    std::vector<unsigned char> code((size_t)n);
    simd::encode_bases(s, code.data(), (size_t)n);

    // Fast path for the common single-window load: the whole seq_ is genuine, so index every position
    // (byte-identical to before). When ensure() created disjoint windows, skip any k-mer that spans a
    // gap so the index matches Java's sparse reference map.
    const bool fullGenuine = genuine_.size() == 1
                             && genuine_[0].first <= loadedStart_ && genuine_[0].second >= loadedEnd();

    auto build = [&](int k, std::unordered_map<uint64_t, int>& m) {
        if (n < k) return;
        const uint64_t lowMask = (1ULL << (3 * (k - 1))) - 1; // keeps all but the top base before shift
        uint64_t enc = 0;
        for (int j = 0; j < k; ++j) enc = (enc << 3) | code[j];  // first window == encodeKmer(s, k)
        auto tryEmplace = [&](int i) {
            int pos = i + loadedStart_;
            if (!fullGenuine && !inGenuine(pos, pos + k - 1)) return;
            auto r = m.emplace(enc, pos);
            if (!r.second) r.first->second = 0;
        };
        tryEmplace(0);
        for (int i = 1; i + k <= n; ++i) {
            enc = ((enc & lowMask) << 3) | code[i + k - 1];     // drop oldest base, append the new one
            tryEmplace(i);
        }
    };
    build(SEED_1, seed17_);
    build(SEED_2, seed12_);
}

void Reference::buildSeedDisjoint() const {
    // Multi-window case: ensure() pulled in disjoint segments for far-off SV/large-indel breakpoints.
    // Index exactly the k-mers that lie fully inside a single genuine window -- byte-identical to the
    // old contiguous-buffer scan filtered by inGenuine(), but reading bases from the disjoint segments
    // via at() instead of materialising the (multi-Mbp) gap between them.
    //
    // Merge the genuine windows into disjoint, position-sorted intervals so every position is visited
    // exactly once: overlapping genuine windows must not double-emplace a k-mer, which would wrongly
    // mark a genuinely-unique k-mer as repeated. Every position inside a merged interval is covered by
    // some genuine window (hence a loaded segment), so at() returns a real base throughout.
    std::vector<std::pair<int,int>> iv = genuine_;
    std::sort(iv.begin(), iv.end());
    std::vector<std::pair<int,int>> merged;
    for (const auto& w : iv) {
        if (!merged.empty() && w.first <= merged.back().second + 1) {
            if (w.second > merged.back().second) merged.back().second = w.second;
        } else {
            merged.push_back(w);
        }
    }
    auto code = [](char ch) -> unsigned {
        return ch=='A'?0u : ch=='C'?1u : ch=='G'?2u : ch=='T'?3u : 4u;  // matches encodeKmer()
    };
    auto build = [&](int k, std::unordered_map<uint64_t, int>& m) {
        const uint64_t lowMask = (1ULL << (3 * (k - 1))) - 1;
        for (const auto& w : merged) {
            int lo = w.first, hi = w.second;
            if (hi - lo + 1 < k) continue;
            uint64_t enc = 0;
            for (int j = 0; j < k; ++j) enc = (enc << 3) | code(at(lo + j));   // == encodeKmer(bases, k)
            auto tryEmplace = [&](int i) {
                if (!inGenuine(i, i + k - 1)) return;   // only k-mers fully inside a single genuine window
                auto r = m.emplace(enc, i);
                if (!r.second) r.first->second = 0;
            };
            tryEmplace(lo);
            for (int i = lo + 1; i + k - 1 <= hi; ++i) {
                enc = ((enc & lowMask) << 3) | code(at(i + k - 1));           // drop oldest, append new
                tryEmplace(i);
            }
        }
    };
    build(SEED_1, seed17_);
    build(SEED_2, seed12_);
}

} // namespace vardict
