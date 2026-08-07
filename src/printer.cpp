#include "printer.hpp"
#include "fisher.hpp"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

namespace vardict {

// Java Utils.getRoundedValueToPrint: integral values print as "0"-pattern (no decimals), otherwise
// DecimalFormat(pattern).format(value) with trailing zeros stripped (may leave a trailing dot,
// exactly as Java does). `decimals` is the number of fraction digits in the pattern.
// Java Utils.roundHalfEven(pattern, value): DecimalFormat(HALF_EVEN).format then parse back.
static double roundHalfEven(int decimals, double value) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, value);
    return std::atof(buf);
}

static std::string getRounded(int decimals, double v) {
    char buf[64];
    if (v == std::floor(v + 0.5)) { // value == Math.round(value): integral
        std::snprintf(buf, sizeof(buf), "%.0f", v);
        return buf;
    }
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, v);
    std::string s(buf);
    size_t last = s.find_last_not_of('0');
    s.erase(last + 1); // strip trailing zeros (replaceAll("0+$",""))
    return s;
}

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

// 38-column fisher-mode line (SimpleOutputVariant.create_simple_variant_38columns): inserts the
// strand-bias Fisher p-value + odds-ratio after QStd, and switches numeric formatting to
// getRoundedValueToPrint (drops trailing .0 / trailing zeros).
static void appendVariantFisher(std::string& out, const Config& cfg, const Region& region, const Variant& v) {
    FisherExact fisher(v.refFwd, v.refRev, v.varFwd, v.varRev);
    std::string pvalue   = getRounded(5, fisher.getPValue());
    std::string oddratio = fisher.getOddRatio();

    // ToVarsBuilder pre-rounds these before either printer: frequency to 4 decimals, mean
    // position/quality/mapping-quality to 1 decimal (roundHalfEven). getRounded then drops the
    // trailing ".0" only when the pre-rounded value is truly integral, matching Java exactly.
    std::string af    = getRounded(4, roundHalfEven(4, v.frequency));
    std::string pmean = getRounded(1, roundHalfEven(1, v.pmean));
    std::string qmean = getRounded(1, roundHalfEven(1, v.qmean));
    std::string mq    = getRounded(1, roundHalfEven(1, v.mapq));
    std::string sn    = getRounded(3, v.qratio);
    // hifreq/nm use fixed DecimalFormat (not stripped), "0" when zero.
    char hbuf[64]; std::string hiaf;
    if (v.hifreq == 0) hiaf = "0"; else { std::snprintf(hbuf, sizeof(hbuf), "%.4f", v.hifreq); hiaf = hbuf; }
    std::string exaf  = getRounded(4, v.extrafreq);
    std::string msi   = getRounded(3, v.msi);
    std::string nm;
    if (v.nm > 0) { char nbuf[64]; std::snprintf(nbuf, sizeof(nbuf), "%.1f", v.nm); nm = nbuf; } else nm = "0";
    std::string dup   = getRounded(2, v.duprate);

    // Size to the variable-length content (see appendVariant); a fixed buffer truncates large complex
    // variants, losing the trailing columns and newline and merging the next variant onto the line.
    size_t cap = 512 + cfg.sample.size() + region.gene.size() + 2 * region.chr.size()
               + v.refallele.size() + v.varallele.size() + v.genotype.size() + v.bias.size()
               + v.leftseq.size() + v.rightseq.size() + v.vartype.size() + v.svInfo.size()
               + pvalue.size() + oddratio.size();
    std::vector<char> buf(cap);
    std::snprintf(buf.data(), cap,
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t"
        "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s:%d-%d\t%s\t%s\t%s\n",
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
        pvalue.c_str(),
        oddratio.c_str(),
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
        v.svInfo.empty() ? "0" : v.svInfo.c_str());  // SV_info
    out += buf.data();
}

void appendVariant(std::string& out, const Config& cfg, const Region& region, const Variant& v) {
    if (cfg.fisher) { appendVariantFisher(out, cfg, region, v); return; }
    // Column order matches SimpleOutputVariant.create_simple_variant_36columns().
    std::string af    = fmt(v.frequency, "%.4f");
    std::string pmean = fmt(v.pmean, "%.1f");
    std::string qmean = fmt(v.qmean, "%.1f");
    // Java stores meanMappingQuality/numberOfMismatches pre-rounded (ToVarsBuilder roundHalfEven "0.0")
    // and then prints "0" when that STORED value is 0 (mapq) / not > 0 (nm). cpp keeps the raw mean, so
    // base the zero test on the rounded value; nonzero values still format the raw mean via %.1f (equal
    // to formatting the rounded value). Otherwise a tiny-nonzero mean that rounds to 0.0 prints "0.0".
    std::string mq    = roundHalfEven(1, v.mapq) == 0 ? "0" : fmt(v.mapq, "%.1f");
    std::string sn    = fmt(v.qratio, "%.3f");
    std::string hiaf  = fmt(v.hifreq, "%.4f");
    // ExtraAF: ToVarsBuilder pre-rounds extraFrequency to 4 decimals (roundHalfEven "0.0000",
    // ToVarsBuilder line 921); the column is then emitted at a fixed 4 decimals ("0" only when the
    // rounded value is exactly 0). A tiny extracnt/depth ratio thus rounds to 0.0 and prints "0"
    // (not "0.0000"), while a genuine value keeps its trailing zeros ("0.2500", not "0.25").
    std::string exaf  = fmt(roundHalfEven(4, v.extrafreq), "%.4f");
    std::string msi   = fmt(v.msi, "%.3f");
    // NM prints with 0.0 only when > 0, else "0".
    std::string nm    = roundHalfEven(1, v.nm) > 0 ? fmt(v.nm, "%.1f") : "0";
    std::string dup   = fmt(v.duprate, "%.1f");

    // The refallele/varallele (and left/right seq, genotype) of a large complex/insertion variant can be
    // hundreds of bases, so a fixed buffer would truncate the line -- dropping its trailing columns AND
    // the newline, which merges the following variant onto this line and loses it. Size the buffer to the
    // actual variable-length content plus generous slack for the bounded numeric/format columns.
    size_t cap = 512 + cfg.sample.size() + region.gene.size() + 2 * region.chr.size()
               + v.refallele.size() + v.varallele.size() + v.genotype.size() + v.bias.size()
               + v.leftseq.size() + v.rightseq.size() + v.vartype.size() + v.svInfo.size();
    std::vector<char> buf(cap);
    std::snprintf(buf.data(), cap,
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
        v.svInfo.empty() ? "0" : v.svInfo.c_str());  // SV_info: empty -> "0"
    out += buf.data();
}

} // namespace vardict
