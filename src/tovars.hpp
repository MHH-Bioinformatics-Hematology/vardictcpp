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
    std::string bias = "0";     // matches VarDictJava Variant.strandBiasFlag default; real variants
                                // overwrite this, so only placeholder variants (e.g. a StrongLOH tumor
                                // block) show it -- where Java prints "0", not "0;0".
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
    std::string svInfo;         // "<splits>-<pairs>-<clusters>" for structural variants, else empty ("0")
    // Somatic-mode only: the variant description string (allele map key, used for sort tie-break and
    // cross-sample matching) and whether it passes isGoodVar (with Java-rounded frequency).
    std::string descriptionString;
    bool good = false;
};

// One position's full candidate variant set for somatic (paired) analysis. Unlike simple mode, no
// variant is dropped here (neither the position-level maxfreq gate nor isGoodVar): every candidate is
// kept, sorted as ToVarsBuilder.sortVariants, and flagged good/bad so SomaticPostProcessModule can
// replicate the leading-good-run + LOH logic. refHicnt/refMeanMapq are the reference-allele stats
// needed by isGoodVar for the compared sample.
struct SomaticPosition {
    int position = 0;
    std::vector<Variant> variants;
    int refHicnt = 0;
    double refMeanMapq = 0;
    // Reference-allele variant (mirrors Vars.referenceVariant) and the position SV string, both needed
    // by the two-sample somatic comparison (LOH / StrongLOH varForPrint, strong-somatic fallback).
    Variant referenceVariant;
    bool hasRef = false;
    std::string sv;
};

// Build called variants for a region from its VariationData (mirrors ToVarsBuilder + the
// SimplePostProcessModule filtering, counting subset). Emits reference-context columns
// (shift3/MSI/flanks) at VarDict defaults where realignment is not yet ported.
std::vector<Variant> callVariants(const Config& cfg, const Region& region,
                                  const VariationData& vd, Reference& ref);

// Somatic-mode counterpart: builds the per-position candidate variants for a region WITHOUT the simple
// maxfreq/isGoodVar dropping (mirrors ToVarsBuilder in somatic/hasBam2 mode, where positions whose
// maxfreq <= -f are kept). Each variant carries its good flag and description string; positions are
// returned sorted ascending. Consumed by the somatic post-processing (somatic.cpp).
std::vector<SomaticPosition> callVariantsSomatic(const Config& cfg, const Region& region,
                                                 const VariationData& vd, Reference& ref);

} // namespace vardict
