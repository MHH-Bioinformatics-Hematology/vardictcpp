#include "printer.hpp"
#include <cstdio>

namespace vardict {

void printHeader(std::FILE* out) {
    std::fprintf(out,
        "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\tDepth\tAltDepth\tRefFwdReads\tRefRevReads\t"
        "AltFwdReads\tAltRevReads\tGenotype\tAF\tBias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\t"
        "HiAF\tExtraAF\tshift3\tMSI\tMSI_NT\tNM\tHiCnt\tHiCov\t5pFlankSeq\t3pFlankSeq\tSeg\t"
        "VarType\tDuprate\tSV_info\n");
}

void printVariant(std::FILE* out, const Config& cfg, const Region& region, const Variant& v) {
    // Column order matches SimpleOutputVariant.create_simple_variant_36columns().
    std::fprintf(out,
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%.4f\t%s\t%.1f\t%d\t%.1f\t%d\t"
        "%.1f\t%.3f\t%.4f\t%.4f\t%d\t%.3f\t%d\t%.1f\t%d\t%d\t%s\t%s\t%s:%d-%d\t%s\t%.4f\t\n",
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
        v.frequency,
        v.bias.c_str(),
        v.pmean,
        v.pstd,
        v.qmean,
        v.qstd,
        v.mapq,
        v.qratio,
        v.hifreq,
        v.extrafreq,
        v.shift3,
        v.msi,
        v.msint,
        v.nm,
        v.hicnt,
        v.hicov,
        v.leftseq.empty() ? "0" : v.leftseq.c_str(),
        v.rightseq.empty() ? "0" : v.rightseq.c_str(),
        region.chr.c_str(), region.start, region.end,
        v.vartype.c_str(),
        v.duprate);
}

} // namespace vardict
