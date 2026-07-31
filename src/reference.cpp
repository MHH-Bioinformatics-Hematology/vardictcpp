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
    if (n > 0) seed17_.reserve((size_t)n);
    const char* s = seq_.data();
    for (int i = 0; i + SEED_1 <= n; ++i) {
        auto r = seed17_.emplace(encodeKmer(s + i, SEED_1), i + loadedStart_);
        if (!r.second) r.first->second = 0;               // second occurrence -> not unique
    }
    for (int i = 0; i + SEED_2 <= n; ++i) {
        auto r = seed12_.emplace(encodeKmer(s + i, SEED_2), i + loadedStart_);
        if (!r.second) r.first->second = 0;
    }
}

} // namespace vardict
