#pragma once
#include <string>
#include <vector>
#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "cigar_parser.hpp"

namespace vardict {

// Output-ready variant at a position (subset of variations/Variant.java + printers/SimpleOutputVariant).
struct Variant {
    std::string refallele;
    std::string varallele;
    int startPosition = 0;
    int endPosition = 0;
    int totalPosCoverage = 0;   // Depth
    int varsCount = 0;          // AltDepth
    int refFwd = 0, refRev = 0; // reference strand reads
    int varFwd = 0, varRev = 0; // variant strand reads
    std::string genotype;
    double frequency = 0;       // AF
    std::string bias = "0;0";
    double pmean = 0;           // PMean
    int    pstd = 0;            // PStd (0/1)
    double qmean = 0;           // QMean
    int    qstd = 0;            // QStd (0/1)
    double mapq = 0;            // MQ
    double qratio = 0;          // Sig_Noise
    double hifreq = 0;          // HiAF
    double extrafreq = 0;       // ExtraAF
    int    shift3 = 0;
    double msi = 0;
    int    msint = 0;
    double nm = 0;              // NM
    int    hicnt = 0;           // HiCnt
    int    hicov = 0;           // HiCov
    std::string leftseq;        // 5' flank
    std::string rightseq;       // 3' flank
    std::string vartype = "SNV";
    double duprate = 0;
};

// Build called variants for a region from its VariationData (mirrors ToVarsBuilder + the
// SimplePostProcessModule filtering, counting subset). Emits reference-context columns
// (shift3/MSI/flanks) at VarDict defaults where realignment is not yet ported.
std::vector<Variant> callVariants(const Config& cfg, const Region& region,
                                  const VariationData& vd, Reference& ref);

} // namespace vardict
