#include "printer.hpp"
#include <cstdio>
#include <string>

namespace vardict {

void printHeader(std::FILE* out) {
    std::fprintf(out,
        "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\tDepth\tAltDepth\tRefFwdReads\tRefRevReads\t"
        "AltFwdReads\tAltRevReads\tGenotype\tAF\tBias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\t"
        "HiAF\tExtraAF\tshift3\tMSI\tMSI_NT\tNM\tHiCnt\tHiCov\t5pFlankSeq\t3pFlankSeq\tSeg\t"
        "VarType\tDuprate\tSV_info\n");
}

// VarDict prints numeric fields as "0" when exactly zero, else with a fixed DecimalFormat pattern.
static std::string fmt(double v, const char* pat) {
    if (v == 0) return "0";
    char buf[64];
    std::snprintf(buf, sizeof(buf), pat, v);
    return buf;
}

void printVariant(std::FILE* out, const Config& cfg, const Region& region, const Variant& v) {
    // Column order matches SimpleOutputVariant.create_simple_variant_36columns().
    std::string af    = fmt(v.frequency, "%.4f");
    std::string pmean = fmt(v.pmean, "%.1f");
    std::string qmean = fmt(v.qmean, "%.1f");
    std::string mq    = fmt(v.mapq, "%.1f");
    std::string sn    = fmt(v.qratio, "%.3f");
    std::string hiaf  = fmt(v.hifreq, "%.4f");
    std::string exaf  = fmt(v.extrafreq, "%.4f");
    std::string msi   = fmt(v.msi, "%.3f");
    // NM prints with 0.0 only when > 0, else "0".
    std::string nm    = v.nm > 0 ? fmt(v.nm, "%.1f") : "0";
    std::string dup   = fmt(v.duprate, "%.1f");

    std::fprintf(out,
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t"
        "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s:%d-%d\t%s\t%s\t%s\n",
        cfg.sample.c_str(),
        region.gene.c_str(),
        region.chr.c_str(),
        v.startPosition,
        v.endPosition,
        v.refallele.c_str(),
        v.varallele.c_str(),
        v.totalPosCoverage,
        v.varsCount,
        v.refFwd, v.refRev,
        v.varFwd, v.varRev,
        v.genotype.c_str(),
        af.c_str(),
        v.bias.c_str(),
        pmean.c_str(),
        v.pstd,
        qmean.c_str(),
        v.qstd,
        mq.c_str(),
        sn.c_str(),
        hiaf.c_str(),
        exaf.c_str(),
        v.shift3,
        msi.c_str(),
        v.msint,
        nm.c_str(),
        v.hicnt,
        v.hicov,
        v.leftseq.empty() ? "0" : v.leftseq.c_str(),
        v.rightseq.empty() ? "0" : v.rightseq.c_str(),
        region.chr.c_str(), region.start, region.end,
        v.vartype.c_str(),
        dup.c_str(),
        "0");  // SV_info: empty -> "0"
}

} // namespace vardict
