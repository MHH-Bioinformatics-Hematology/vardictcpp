#pragma once
#include <string>
#include <unordered_map>
#include <map>
#include <cstdint>

namespace vardict {

// Mirrors VarDictJava variations/Variation.java. Per-allele accumulator at a position.
struct Variation {
    int   varsCount = 0;            // $cnt
    int   varsCountOnForward = 0;   // $dirPlus
    int   varsCountOnReverse = 0;   // $dirMinus
    double meanPosition = 0;        // $pmean  (sum of read positions)
    double meanQuality = 0;         // $qmean  (sum base qualities)
    double meanMappingQuality = 0;  // $Qmean  (sum mapping qualities)
    double numberOfMismatches = 0;  // $nm     (sum NM)
    int   lowQualityReadsCount = 0; // $locnt
    int   highQualityReadsCount = 0;// $hicnt
    bool  pstd = false;             // seen at >=2 distinct read positions
    bool  qstd = false;             // seen with >=2 distinct qualities
    int   pp = 0;                   // previous read position (for pstd)
    double pq = 0;                  // previous base quality (for qstd)
    int   extracnt = 0;             // indel realignment adjustment (not yet populated)

    inline void incDir(bool reverse) { reverse ? ++varsCountOnReverse : ++varsCountOnForward; }
};

// A position's allele map. std::map keeps a deterministic (ascending key) order for reproducible
// output; VarDict uses insertion-ordered LinkedHashMap and re-sorts before printing, so the
// container choice does not change the emitted set, only intra-position tie handling.
using VarMap = std::map<std::string, Variation>;

} // namespace vardict
