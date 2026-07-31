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

void Reference::load(const std::string& chr, int start, int end, int pad) {
    int s = start - pad; if (s < 1) s = 1;
    int e = end + pad;
    int len = 0;
    // faidx_fetch_seq is 0-based, end-inclusive.
    char* raw = faidx_fetch_seq(fai_, chr.c_str(), s - 1, e - 1, &len);
    if (!raw || len < 0) {
        // Region outside contig / not found: leave empty so at() returns 'N'.
        seq_.clear();
        loadedStart_ = s;
        seedBuilt_ = false;
        if (raw) free(raw);
        return;
    }
    seq_.assign(raw, len);
    // Upper-case the whole loaded window (up to a full chromosome) with SIMD; identical result to the
    // per-byte std::toupper, just 16 bases at a time. Portable across SSE2 / NEON / scalar.
    if (!seq_.empty()) simd::toupper_ascii(&seq_[0], seq_.size());
    loadedStart_ = s;
    loadedChr_ = chr;
    seedBuilt_ = false;
    free(raw);
}

void Reference::ensure(int start, int end) {
    if (loadedChr_.empty()) return;
    if (start >= loadedStart_ && end <= loadedEnd()) return; // already covered
    int s = std::min(start, loadedStart_); if (s < 1) s = 1;
    int e = std::max(end, loadedEnd());
    load(loadedChr_, s, e, 0);
}

void Reference::buildSeed() const {
    seed17_.clear(); seed12_.clear();
    seedBuilt_ = true;
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

    auto build = [&](int k, std::unordered_map<uint64_t, int>& m) {
        if (n < k) return;
        const uint64_t lowMask = (1ULL << (3 * (k - 1))) - 1; // keeps all but the top base before shift
        uint64_t enc = 0;
        for (int j = 0; j < k; ++j) enc = (enc << 3) | code[j];  // first window == encodeKmer(s, k)
        {
            auto r = m.emplace(enc, 0 + loadedStart_);
            if (!r.second) r.first->second = 0;
        }
        for (int i = 1; i + k <= n; ++i) {
            enc = ((enc & lowMask) << 3) | code[i + k - 1];     // drop oldest base, append the new one
            auto r = m.emplace(enc, i + loadedStart_);
            if (!r.second) r.first->second = 0;
        }
    };
    build(SEED_1, seed17_);
    build(SEED_2, seed12_);
}

} // namespace vardict
