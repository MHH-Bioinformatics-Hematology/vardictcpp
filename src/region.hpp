#pragma once
#include <string>
#include <vector>
#include <algorithm>

namespace vardict {

struct Region {
    std::string chr;
    int start = 0;
    int end = 0;
    std::string gene;
    // Amplicon mode only (data/Region.java insertStart/insertEnd = BED thickStart/thickEnd). The insert
    // is the amplicon's on-target interval; positions outside it are primer/flank and not reported.
    int insertStart = 0;
    int insertEnd = 0;
};

// Split regions longer than chunkSize into consecutive windows (mirrors the Java --chunk change),
// bounding per-region memory. chunkSize <= 0 leaves regions unchanged.
inline std::vector<Region> splitLongRegions(const std::vector<Region>& in, int chunkSize) {
    if (chunkSize <= 0) return in;
    std::vector<Region> out;
    for (const auto& r : in) {
        if (r.end - r.start + 1 <= chunkSize) { out.push_back(r); continue; }
        for (int s = r.start; s <= r.end; s += chunkSize) {
            Region w = r;
            w.start = s;
            w.end = std::min(s + chunkSize - 1, r.end);
            out.push_back(w);
        }
    }
    return out;
}

} // namespace vardict
