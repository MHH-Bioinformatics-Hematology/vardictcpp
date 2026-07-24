#pragma once
// Ports variations/Variation.java, Sclip.java, Mate.java.
#include <string>
#include <map>
#include <vector>
#include <memory>

namespace vardict {

// variations/Variation.java -- per-allele accumulator at a position.
struct Variation {
    int   varsCount = 0;
    int   varsCountOnForward = 0;
    int   varsCountOnReverse = 0;
    double meanPosition = 0;
    double meanQuality = 0;
    double meanMappingQuality = 0;
    double numberOfMismatches = 0;
    int   lowQualityReadsCount = 0;
    int   highQualityReadsCount = 0;
    bool  pstd = false;
    bool  qstd = false;
    int   pp = 0;
    double pq = 0;
    int   extracnt = 0;

    inline void incDir(bool reverse) { reverse ? ++varsCountOnReverse : ++varsCountOnForward; }
    inline void decDir(bool reverse) { reverse ? --varsCountOnReverse : --varsCountOnForward; }
    inline int  getDir(bool reverse) const { return reverse ? varsCountOnReverse : varsCountOnForward; }
    inline void addDir(bool reverse, int add) { if (reverse) varsCountOnReverse += add; else varsCountOnForward += add; }
    inline void subDir(bool reverse, int sub) { if (reverse) varsCountOnReverse -= sub; else varsCountOnForward -= sub; }
};

// A position's allele map (VarDict's VariationMap without the SV marker, which is tracked separately
// when SV support lands). std::map keeps a deterministic ordered iteration.
using VarMap = std::map<std::string, Variation>;

// variations/Mate.java
struct Mate {
    int mateStart_ms=0, mateEnd_me=0, mateLength_mlen=0, start_s=0, end_e=0;
    double pmean_rp=0, qmean_q=0, Qmean_Q=0, nm=0;
};

// collection/VariationMap.SV: per-position structural-variant marker (pairs/splits/clusters) that
// ToVarsBuilder joins as "<splits>-<pairs>-<clusters>" into the SV_info output column.
struct SVInfo {
    std::string type;
    int pairs=0, splits=0, clusters=0;
};

// variations/Sclip.java (extends Variation)
struct Sclip : Variation {
    std::map<int, std::map<char,int>> nt;                          // pos -> base -> count
    std::map<int, std::map<char, std::shared_ptr<Variation>>> seq; // pos -> base -> Variation
    std::string sequence;
    bool used = false;
    int start=0, end=0, mstart=0, mend=0, mlen=0, disc=0, softp=0;
    std::map<int,int> soft;
    std::vector<Mate> mates;
};

} // namespace vardict
