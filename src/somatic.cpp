#include "somatic.hpp"
#include <cstdio>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

namespace vardict {

// %.4f is round-half-to-even (IEEE default), matching Java Utils.roundHalfEven / createVariant frequency.
static double round4(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.4f", x);
    return std::atof(buf);
}

// Round to n decimals with round-half-to-even (matches Java Utils.roundHalfEven). ToVarsBuilder
// stores pmean/qmean/mapq/nm at 1 dp and frequency/hifreq/extrafreq at 4 dp, so combineAnalysis reads
// those ALREADY-ROUNDED values; cpp stores them raw and must round at read time to match byte-for-byte.
static double roundNloc(double x, int n) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.*f", n, x);
    return std::atof(buf);
}

// VarDict prints a numeric field as "0" when exactly zero, else with a fixed DecimalFormat pattern.
static std::string fmt(double v, const char* pat) {
    if (v == 0) return "0";
    char buf[64];
    std::snprintf(buf, sizeof(buf), pat, v);
    return buf;
}

// Variant.isNoise as a pure predicate (mirrors somatic determinateType's use).
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

// Variant.adjComplex: trim the common leading/trailing bases of a Complex variant's ref/alt alleles and
// shift start/end accordingly. Ported for the somatic COMPLEX print path.
static void adjComplex(Variant& v) {
    std::string ref = v.refallele, var = v.varallele;
    if (!var.empty() && var[0] == '<') return; // structural
    int n = 0;
    while (n < (int)ref.size() - 1 && n < (int)var.size() - 1 && ref[n] == var[n]) n++;
    if (n > 0) {
        v.startPosition += n;
        ref = ref.substr(n);
        var = var.substr(n);
        v.refallele = ref; v.varallele = var;
    }
    n = 0;
    while ((int)ref.size() - 1 - n > 0 && (int)var.size() - 1 - n > 0
           && ref[ref.size() - 1 - n] == var[var.size() - 1 - n]) n++;
    if (n > 0) {
        v.endPosition -= n;
        v.refallele = ref.substr(0, ref.size() - n);
        v.varallele = var.substr(0, var.size() - n);
    }
}

// SomaticPostProcessModule.determinateType: 5-way classification of a compared variant.
static std::string determinateType(const Config& cfg, const Variant& standardVariant, const Variant& variantToCompare) {
    double sfreq = standardVariant.frequency;
    double vfreq = variantToCompare.frequency;
    int vcov = variantToCompare.varsCount;
    std::string type;
    if (variantToCompare.good) { // isGoodVar
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

// VarDict strand-bias flag (variations/VariationUtils.strandBias), reimplemented locally for the
// combineAnalysis back-subtraction bias string. minBiasReads==0 is treated as 2 (matching tovars call sites).
static int strandBiasLocal(int fwd, int rev, const Config& cfg) {
    int minb = cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads;
    if (fwd + rev <= 12) return (fwd > 0 && rev > 0) ? 2 : 0;
    double tot = fwd + rev;
    bool ok = (fwd / tot >= cfg.bias) && (rev / tot >= cfg.bias) && fwd >= minb && rev >= minb;
    return ok ? 2 : 1;
}

// jregex MINUS_NUM_NUM "-\d\d" with .find() (search): a '-' immediately followed by two digits anywhere.
static bool minusNumNum(const std::string& s) {
    for (size_t i = 0; i + 2 < s.size(); ++i)
        if (s[i] == '-' && std::isdigit((unsigned char)s[i + 1]) && std::isdigit((unsigned char)s[i + 2]))
            return true;
    return false;
}

// SomaticPostProcessModule.combineAnalysis: re-run the pipeline over a widened region on the MERGED
// bam1+bam2 read stream (via the `combine` callback), find the merged variant at (position, desc), and
// if the merged coverage exceeds variant1's by >= minr, back-subtract to fill variant2 (the placeholder
// for the OTHER sample) and reclassify as Germline. Returns "FALSE" (drop), "Germline" (reclassified),
// or "" (unchanged). Mirrors the Java field arithmetic and clamps exactly.
static std::string combineAnalysis(const Config& cfg, const CombineFn& combine,
                                   const Variant& variant1, Variant& variant2,
                                   int position, const std::string& desc) {
    // Don't do it for structural variants.
    if (variant1.endPosition - variant1.startPosition > cfg.SVMINLEN) return "";
    if (!combine) return ""; // callback disabled -> no refinement
    Variant vref;
    if (!combine(variant1.startPosition, variant1.endPosition, position, desc, vref)) {
        return "FALSE"; // no matching merged variant (Java: getVarMaybe == null)
    }
    if (vref.varsCount - variant1.varsCount >= cfg.minReads) {
        variant2.totalPosCoverage = std::max(0, vref.totalPosCoverage - variant1.totalPosCoverage);
        variant2.varsCount        = std::max(0, vref.varsCount        - variant1.varsCount);
        variant2.refFwd           = std::max(0, vref.refFwd           - variant1.refFwd);
        variant2.refRev           = std::max(0, vref.refRev           - variant1.refRev);
        variant2.varFwd           = std::max(0, vref.varFwd           - variant1.varFwd);
        variant2.varRev           = std::max(0, vref.varRev           - variant1.varRev);
        if (variant2.varsCount != 0) {
            double pc = variant2.varsCount;
            // Java reads the stored (rounded) mean fields: pmean/qmean/mapq/nm at 1 dp, hifreq/extrafreq
            // at 4 dp. Round the raw cpp values here so the back-subtraction matches byte-for-byte.
            double v1p = roundNloc(variant1.pmean, 1),     vrp = roundNloc(vref.pmean, 1);
            double v1q = roundNloc(variant1.qmean, 1),     vrq = roundNloc(vref.qmean, 1);
            double v1m = roundNloc(variant1.mapq, 1),      vrm = roundNloc(vref.mapq, 1);
            double v1h = roundNloc(variant1.hifreq, 4),    vrh = roundNloc(vref.hifreq, 4);
            double v1e = roundNloc(variant1.extrafreq, 4), vre = roundNloc(vref.extrafreq, 4);
            double v1n = roundNloc(variant1.nm, 1),        vrn = roundNloc(vref.nm, 1);
            variant2.pmean     = (vrp * vref.varsCount - v1p * variant1.varsCount) / pc;
            variant2.qmean     = (vrq * vref.varsCount - v1q * variant1.varsCount) / pc;
            variant2.mapq      = (vrm * vref.varsCount - v1m * variant1.varsCount) / pc;
            variant2.hifreq    = (vrh * vref.varsCount - v1h * variant1.varsCount) / pc;
            variant2.extrafreq = (vre * vref.varsCount - v1e * variant1.varsCount) / pc;
            variant2.nm        = (vrn * vref.varsCount - v1n * variant1.varsCount) / pc;
        } else {
            variant2.pmean = 0; variant2.qmean = 0; variant2.mapq = 0;
            variant2.hifreq = 0; variant2.extrafreq = 0; variant2.nm = 0;
        }
        variant2.pstd = 1; // isAtLeastAt2Positions = true
        variant2.qstd = 1; // hasAtLeast2DiffQualities = true
        if (variant2.totalPosCoverage <= 0) return "FALSE";
        variant2.frequency = variant2.varsCount / (double)variant2.totalPosCoverage;
        variant2.qratio = variant1.qratio; // Can't back-calculate; inherits variant1's ratio
        variant2.genotype = vref.genotype;
        variant2.bias = std::to_string(strandBiasLocal(variant2.refFwd, variant2.refRev, cfg)) + ";" +
                        std::to_string(strandBiasLocal(variant2.varFwd, variant2.varRev, cfg));
        return "Germline";
    } else if (vref.varsCount < variant1.varsCount - 2) {
        return "FALSE";
    }
    return "";
}

// The 18-field per-sample block (Depth..NM). A null slot prints 18 zeros (matching a null Variant).
static std::string block(const Variant* v) {
    if (!v) return "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
    char buf[1024];
    std::string geno = v->genotype.empty() ? "0" : v->genotype;
    std::string bias = v->bias.empty() ? "0" : v->bias;
    std::snprintf(buf, sizeof(buf),
        "%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s",
        v->totalPosCoverage, v->varsCount, v->refFwd, v->refRev, v->varFwd, v->varRev,
        geno.c_str(), fmt(v->frequency, "%.4f").c_str(), bias.c_str(),
        fmt(v->pmean, "%.1f").c_str(), v->pstd, fmt(v->qmean, "%.1f").c_str(), v->qstd,
        fmt(v->mapq, "%.1f").c_str(), fmt(v->qratio, "%.3f").c_str(),
        fmt(v->hifreq, "%.4f").c_str(), fmt(v->extrafreq, "%.4f").c_str(),
        v->nm > 0 ? fmt(v->nm, "%.1f").c_str() : "0");
    return buf;
}

// SomaticOutputVariant: begin -> pos/ref/alt/vartype, end -> shift3/msi/flanks, tumor -> var1 block,
// normal -> var2 block. sv1/sv2 -> the two SV_info columns; label -> VarLabel (the somatic type).
static void printSomatic(std::string& out, const std::string& sample, const Region& region,
                         const Variant* begin, const Variant* end, const Variant* tumor, const Variant* normal,
                         const std::string& sv1, const std::string& sv2, const std::string& label) {
    if (!begin) return;
    std::string leftS  = (end && !end->leftseq.empty())  ? end->leftseq  : "0";
    std::string rightS = (end && !end->rightseq.empty()) ? end->rightseq : "0";
    char head[512];
    std::snprintf(head, sizeof(head), "%s\t%s\t%s\t%d\t%d\t%s\t%s\t",
        sample.c_str(), region.gene.c_str(), region.chr.c_str(),
        begin->startPosition, begin->endPosition, begin->refallele.c_str(), begin->varallele.c_str());
    char tail[512];
    std::snprintf(tail, sizeof(tail), "\t%d\t%s\t%d\t%s\t%s\t%s:%d-%d\t%s\t%s\t%s\t%s\t%s\t%s\n",
        end ? end->shift3 : 0, end ? fmt(end->msi, "%.3f").c_str() : "0", end ? end->msint : 0,
        leftS.c_str(), rightS.c_str(),
        region.chr.c_str(), region.start, region.end,
        label.c_str(), begin->vartype.c_str(),
        tumor ? fmt(tumor->duprate, "%.1f").c_str() : "0", sv1.empty() ? "0" : sv1.c_str(),
        normal ? fmt(normal->duprate, "%.1f").c_str() : "0", sv2.empty() ? "0" : sv2.c_str());
    out += head;
    out += block(tumor);
    out += "\t";
    out += block(normal);
    out += tail;
}

// Lookups mirroring getVarMaybe(vars, varn, nt) and getVarMaybe(vars, var, 0).
static const Variant* getVarByDesc(const SomaticPosition* sp, const std::string& nt) {
    if (!sp) return nullptr;
    for (const Variant& v : sp->variants) if (v.descriptionString == nt) return &v;
    return nullptr;
}
static const Variant* getTopVar(const SomaticPosition* sp) {
    if (!sp || sp->variants.empty()) return nullptr;
    return &sp->variants[0];
}

// SomaticPostProcessModule.callingForOneSample: variants exist in only one sample.
static void callingForOneSample(std::string& out, const Config& cfg, const std::string& sample,
                                const Region& region, const SomaticPosition* variants,
                                bool isFirstCover, const std::string& varLabel) {
    if (!variants || variants->variants.empty()) return;
    for (const Variant& variant0 : variants->variants) {
        if (variant0.refallele == variant0.varallele) continue;
        if (!variant0.good) continue; // isGoodVar
        Variant variant = variant0;
        if (variant.vartype == "Complex") adjComplex(variant);
        if (isFirstCover) {
            printSomatic(out, sample, region, &variant, &variant, nullptr, &variant, "", variants->sv, varLabel);
        } else {
            printSomatic(out, sample, region, &variant, &variant, &variant, nullptr, variants->sv, "", varLabel);
        }
    }
}

// SomaticPostProcessModule.printVariationsFromFirstSample.
static void printVariationsFromFirstSample(std::string& out, const Config& cfg, const std::string& sample,
                                           const Region& region, int position,
                                           const SomaticPosition* v1, const SomaticPosition* v2,
                                           const CombineFn& combine) {
    const std::string& sv1 = v1->sv;
    const std::string& sv2 = v2->sv;
    size_t n = 0;
    while (n < v1->variants.size() && v1->variants[n].good) {
        const Variant& vref0 = v1->variants[n];
        if (vref0.refallele == vref0.varallele) { ++n; continue; }
        std::string nt = vref0.descriptionString;
        Variant vref = vref0;
        if (vref.vartype == "Complex") adjComplex(vref);
        const Variant* v2nt = getVarByDesc(v2, nt);
        if (v2nt != nullptr) {
            std::string type = determinateType(cfg, vref, *v2nt);
            printSomatic(out, sample, region, &vref, v2nt, &vref, v2nt, sv1, sv2, type);
        } else { // sample 1 only, should be strong somatic
            Variant varForPrint;
            const Variant* varForPrintPtr = nullptr;
            if (!v2->variants.empty()) {
                const Variant* v2r = getTopVar(v2);
                varForPrint.totalPosCoverage = (v2r && v2r->totalPosCoverage) ? v2r->totalPosCoverage : 0;
                varForPrint.refFwd = (v2r && v2r->refFwd) ? v2r->refFwd : 0;
                varForPrint.refRev = (v2r && v2r->refRev) ? v2r->refRev : 0;
                varForPrintPtr = &varForPrint;
            } else if (v2->hasRef) {
                varForPrintPtr = &v2->referenceVariant;
            } else {
                varForPrintPtr = nullptr;
            }
            std::string type = "StrongSomatic";
            Variant v2nt2; // Java: new Variant() placeholder, filled by combineAnalysis on Germline
            if (vref.vartype != "SNV" && (nt.size() > 10 || minusNumNum(nt))) {
                if (vref.varsCount < cfg.minReads + 3 && nt.find('<') == std::string::npos) {
                    std::string newtype = combineAnalysis(cfg, combine, vref, v2nt2, position, nt);
                    if (newtype == "FALSE") { ++n; continue; }
                    if (!newtype.empty()) type = newtype;
                }
            }
            if (type == "StrongSomatic") {
                printSomatic(out, sample, region, &vref, &vref, &vref, varForPrintPtr, sv1, sv2, "StrongSomatic");
            } else {
                printSomatic(out, sample, region, &vref, &vref, &vref, &v2nt2, sv1, sv2, type);
            }
        }
        ++n;
    }
    if (n == 0) {
        if (v2->variants.empty()) return;
        for (const Variant& v2var0 : v2->variants) {
            if (!v2var0.good) continue;
            std::string nt = v2var0.descriptionString;
            const Variant* v1nt = getVarByDesc(v1, nt);
            if (v1nt != nullptr) {
                if (v1nt->refallele == v1nt->varallele) continue;
                Variant v1ntc = *v1nt;
                Variant v2varc = v2var0;
                if (v2varc.vartype == "Complex") adjComplex(v1ntc);
                std::string type = v1ntc.frequency < cfg.lofreq ? "LikelyLOH" : "Germline";
                printSomatic(out, sample, region, &v1ntc, &v2varc, &v1ntc, &v2varc, sv1, sv2, type);
            } else {
                if (v2var0.refallele == v2var0.varallele) continue;
                const Variant* v1var = getTopVar(v1);
                int tcov = (v1var && v1var->totalPosCoverage) ? v1var->totalPosCoverage : 0;
                const Variant* v1ref = v1->hasRef ? &v1->referenceVariant : nullptr;
                int fwd = v1ref ? v1ref->varFwd : 0;
                int rev = v1ref ? v1ref->varRev : 0;
                std::string genotype = v1var ? v1var->genotype
                                             : (v1ref ? v1ref->descriptionString + "/" + v1ref->descriptionString : "N/N");
                Variant v2varc = v2var0;
                if (v2varc.vartype == "Complex") adjComplex(v2varc);
                Variant varForPrint;
                varForPrint.totalPosCoverage = tcov;
                varForPrint.refFwd = fwd;
                varForPrint.refRev = rev;
                varForPrint.genotype = genotype;
                printSomatic(out, sample, region, &v2varc, &v2varc, &varForPrint, &v2varc, "", sv2, "StrongLOH");
            }
        }
    }
}

// SomaticPostProcessModule.printVariationsFromSecondSample (sample 1 has only reference).
static void printVariationsFromSecondSample(std::string& out, const Config& cfg, const std::string& sample,
                                            const Region& region, int position,
                                            const SomaticPosition* v1, const SomaticPosition* v2,
                                            const CombineFn& combine) {
    const std::string& sv2 = v2->sv;
    for (const Variant& v2var0 : v2->variants) {
        if (v2var0.refallele == v2var0.varallele) continue;
        if (!v2var0.good) continue;
        const std::string& descriptionString = v2var0.descriptionString;
        std::string type = "StrongLOH";
        Variant v1nt; // Java: v1.varDescriptionStringToVariants.computeIfAbsent -> new Variant() (posCov 0)
        std::string newType;
        if (v2var0.varsCount < cfg.minReads + 3 && descriptionString.find('<') == std::string::npos
                && (descriptionString.size() > 10 || minusNumNum(descriptionString))) {
            newType = combineAnalysis(cfg, combine, v2var0, v1nt, position, descriptionString);
            if (newType == "FALSE") continue;
        }
        const Variant* varForPrint;
        if (!newType.empty()) {
            type = newType;
            varForPrint = &v1nt;
        } else {
            varForPrint = v1->hasRef ? &v1->referenceVariant : nullptr;
        }
        Variant v2varc = v2var0;
        if (v2varc.vartype == "Complex") adjComplex(v2varc);
        printSomatic(out, sample, region, &v2varc, &v2varc, varForPrint, &v2varc, "", sv2, type);
    }
}

static void callingForBothSamples(std::string& out, const Config& cfg, const std::string& sample,
                                  const Region& region, int position,
                                  const SomaticPosition* v1, const SomaticPosition* v2,
                                  const CombineFn& combine) {
    if (v1->variants.empty() && v2->variants.empty()) return;
    if (!v1->variants.empty()) {
        printVariationsFromFirstSample(out, cfg, sample, region, position, v1, v2, combine);
    } else if (!v2->variants.empty()) {
        printVariationsFromSecondSample(out, cfg, sample, region, position, v1, v2, combine);
    }
}

void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& tumor,
                         const std::vector<SomaticPosition>& normal,
                         const CombineFn& combine) {
    std::string sample = cfg.sample;
    if (!cfg.sample2.empty()) sample += "|" + cfg.sample2;

    std::map<int, const SomaticPosition*> m1, m2;
    for (const SomaticPosition& sp : tumor) m1[sp.position] = &sp;
    for (const SomaticPosition& sp : normal) m2[sp.position] = &sp;

    std::set<int> positions;
    for (const auto& kv : m1) positions.insert(kv.first);
    for (const auto& kv : m2) positions.insert(kv.first);

    for (int position : positions) {
        if (position < region.start || position > region.end) continue;
        auto i1 = m1.find(position); const SomaticPosition* v1 = i1 == m1.end() ? nullptr : i1->second;
        auto i2 = m2.find(position); const SomaticPosition* v2 = i2 == m2.end() ? nullptr : i2->second;
        if (!v1 && !v2) continue;
        if (!v1) {
            callingForOneSample(out, cfg, sample, region, v2, true, "Deletion");
        } else if (!v2) {
            callingForOneSample(out, cfg, sample, region, v1, false, "SampleSpecific");
        } else {
            callingForBothSamples(out, cfg, sample, region, position, v1, v2, combine);
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
