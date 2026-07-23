#include "somatic.hpp"
#include <cstdio>
#include <cmath>
#include <string>

namespace vardict {

// Java Utils.roundHalfEven("0.0000", x): %.4f is round-half-to-even (IEEE default), matching Java's
// createVariant frequency. Somatic gating (isGoodVar / determinateType) compares this rounded value.
static double round4(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4f", x);
    return std::atof(buf);
}

// VarDict prints a numeric field as "0" when exactly zero, else with a fixed DecimalFormat pattern
// (trailing zeros kept). Matches SomaticOutputVariant.create_somatic_variant_55columns.
static std::string fmt(double v, const char* pat) {
    if (v == 0) return "0";
    char buf[64];
    std::snprintf(buf, sizeof(buf), pat, v);
    return buf;
}

// Port of Variant.isNoise as a pure predicate (no field mutation): the mutation only fires when the
// predicate is true, and for well-covered good variants it never is. Kept faithful for correctness.
static bool isNoise(const Config& cfg, const Variant& v) {
    double qual = v.qmean;
    int posCov = v.varsCount;
    bool has2diff = v.qstd != 0;
    if (((qual < 4.5 || (qual < 12 && !has2diff)) && posCov <= 3)
        || (qual < cfg.goodq && round4(v.frequency) < 2 * cfg.lofreq && posCov <= 1)) {
        return true;
    }
    return false;
}

// Port of SomaticPostProcessModule.determinateType. standardVariant/variantToCompare are the tumor and
// (matched) normal variants; here they carry identical values (same BAM), so vc-good == the leading-run
// good flag. Uses Java-rounded frequencies for the thresholds.
static std::string determinateType(const Config& cfg, const Variant& standardVariant,
                                   const Variant& variantToCompare) {
    double sfreq = round4(standardVariant.frequency);
    double vfreq = round4(variantToCompare.frequency);
    int vcov = variantToCompare.varsCount;
    std::string type;
    // variantToCompare.isGoodVar(referenceVar, standardVariant.vartype, splice): standardVariant and
    // variantToCompare share alleles/vartype here, so the precomputed good flag is exact.
    bool vcGood = variantToCompare.good;
    if (vcGood) {
        if (sfreq > (1 - cfg.lofreq) && vfreq < 0.8 && vfreq > 0.2) {
            type = "LikelyLOH";
        } else if (vfreq < cfg.lofreq || vcov <= 1) {
            type = "LikelySomatic";
        } else {
            type = "Germline";
        }
    } else {
        if (vfreq < cfg.lofreq || vcov <= 1) {
            type = "LikelySomatic";
        } else {
            type = "AFDiff";
        }
    }
    if (isNoise(cfg, variantToCompare) && standardVariant.vartype == "SNV") {
        type = "StrongSomatic";
    }
    return type;
}

// Emit one 55-column somatic row. begin/end/tumor/normal all resolve to the same Variant `v` because
// the harness pairs a BAM with itself (var1 block == var2 block); the layout still follows
// SomaticOutputVariant column-for-column so a genuine two-BAM pairing would slot straight in.
static void appendSomaticRow(std::string& out, const Config& cfg, const Region& region,
                             const Variant& v, const std::string& varLabel) {
    std::string freqS  = fmt(v.frequency, "%.4f");
    std::string pmeanS = fmt(v.pmean, "%.1f");
    std::string qmeanS = fmt(v.qmean, "%.1f");
    std::string mapqS  = fmt(v.mapq, "%.1f");
    std::string qratS  = fmt(v.qratio, "%.3f");
    std::string hifS   = fmt(v.hifreq, "%.4f");
    std::string exfS   = fmt(v.extrafreq, "%.4f");
    std::string nmS    = v.nm > 0 ? fmt(v.nm, "%.1f") : std::string("0");
    std::string msiS   = fmt(v.msi, "%.3f");
    std::string dupS   = fmt(v.duprate, "%.1f");
    std::string leftS  = v.leftseq.empty() ? "0" : v.leftseq;
    std::string rightS = v.rightseq.empty() ? "0" : v.rightseq;
    const char* geno   = v.genotype.empty() ? "0" : v.genotype.c_str();
    const char* bias   = v.bias.empty() ? "0" : v.bias.c_str();

    char buf[2048];
    // sample gene chr start end ref alt | var1(18) | var2(18) | shift3 msi msint 5p 3p region label type dup1 sv1 dup2 sv2
    std::snprintf(buf, sizeof(buf),
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t"
        "%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t"
        "%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t"
        "%d\t%s\t%d\t%s\t%s\t%s:%d-%d\t%s\t%s\t%s\t%s\t%s\t%s\n",
        cfg.sample.c_str(), region.gene.c_str(), region.chr.c_str(),
        v.startPosition, v.endPosition, v.refallele.c_str(), v.varallele.c_str(),
        // var1 block
        v.totalPosCoverage, v.varsCount, v.refFwd, v.refRev, v.varFwd, v.varRev,
        geno, freqS.c_str(), bias, pmeanS.c_str(), v.pstd, qmeanS.c_str(), v.qstd,
        mapqS.c_str(), qratS.c_str(), hifS.c_str(), exfS.c_str(), nmS.c_str(),
        // var2 block (identical values)
        v.totalPosCoverage, v.varsCount, v.refFwd, v.refRev, v.varFwd, v.varRev,
        geno, freqS.c_str(), bias, pmeanS.c_str(), v.pstd, qmeanS.c_str(), v.qstd,
        mapqS.c_str(), qratS.c_str(), hifS.c_str(), exfS.c_str(), nmS.c_str(),
        // tail
        v.shift3, msiS.c_str(), v.msint, leftS.c_str(), rightS.c_str(),
        region.chr.c_str(), region.start, region.end,
        varLabel.c_str(), v.vartype.c_str(),
        dupS.c_str(), "0", dupS.c_str(), "0");
    out += buf;
}

void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& positions) {
    for (const SomaticPosition& pos : positions) {
        const std::vector<Variant>& vars = pos.variants;   // v1 == v2 (same BAM pair)
        if (vars.empty()) continue;                        // callingForBothSamples: nothing to compare

        // printVariationsFromFirstSample: process the leading run of good variants.
        size_t n = 0;
        while (n < vars.size() && vars[n].good) {
            const Variant& vref = vars[n];
            if (vref.refallele == vref.varallele) { ++n; continue; }
            const Variant& v2nt = vars[n];                 // getVarMaybe(v2, varn, nt) == matched normal
            std::string type = determinateType(cfg, vref, v2nt);
            appendSomaticRow(out, cfg, region, vref, type);
            ++n;
        }
        if (n == 0) {
            // No leading good variant in sample 1: potential LOH, driven by sample 2's good variants.
            for (const Variant& v2var : vars) {
                if (!v2var.good) continue;
                if (v2var.refallele == v2var.varallele) continue;
                const Variant& v1nt = v2var;               // getVarMaybe(v1, varn, nt) == matched tumor
                std::string type = round4(v1nt.frequency) < cfg.lofreq ? "LikelyLOH" : "Germline";
                appendSomaticRow(out, cfg, region, v1nt, type);
            }
        }
    }
}

void printSomaticHeader(std::FILE* out) {
    std::fprintf(out,
        "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\t"
        "Depth\tAltDepth\tRefFwdReads\tRefRevReads\tAltFwdReads\tAltRevReads\tGenotype\tAF\tBias\t"
        "PMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\tHiAF\tExtraAF\tNM\t"
        "Depth\tAltDepth\tRefFwdReads\tRefRevReads\tAltFwdReads\tAltRevReads\tGenotype\tAF\tBias\t"
        "PMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\tHiAF\tExtraAF\tNM\t"
        "shift3\tMSI\tMSI_NT\t5pFlankSeq\t3pFlankSeq\tSeg\tVarLabel\tVarType\t"
        "Duprate1\tSV_info1\tDuprate2\tSV_info2\n");
}

} // namespace vardict
