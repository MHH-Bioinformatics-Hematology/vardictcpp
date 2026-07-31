#include "somatic.hpp"
#include <cstdio>
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
                                           const Region& region, const SomaticPosition* v1, const SomaticPosition* v2) {
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
        } else { // sample 1 only, strong somatic (combineAnalysis deferred)
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
            printSomatic(out, sample, region, &vref, &vref, &vref, varForPrintPtr, sv1, sv2, "StrongSomatic");
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
                                            const Region& region, const SomaticPosition* v1, const SomaticPosition* v2) {
    const std::string& sv2 = v2->sv;
    for (const Variant& v2var0 : v2->variants) {
        if (v2var0.refallele == v2var0.varallele) continue;
        if (!v2var0.good) continue;
        std::string type = "StrongLOH"; // combineAnalysis deferred
        const Variant* v1ref = v1->hasRef ? &v1->referenceVariant : nullptr;
        const Variant* varForPrint = v1ref;
        Variant v2varc = v2var0;
        if (v2varc.vartype == "Complex") adjComplex(v2varc);
        printSomatic(out, sample, region, &v2varc, &v2varc, varForPrint, &v2varc, "", sv2, type);
    }
}

static void callingForBothSamples(std::string& out, const Config& cfg, const std::string& sample,
                                  const Region& region, const SomaticPosition* v1, const SomaticPosition* v2) {
    if (v1->variants.empty() && v2->variants.empty()) return;
    if (!v1->variants.empty()) {
        printVariationsFromFirstSample(out, cfg, sample, region, v1, v2);
    } else if (!v2->variants.empty()) {
        printVariationsFromSecondSample(out, cfg, sample, region, v1, v2);
    }
}

void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& tumor,
                         const std::vector<SomaticPosition>& normal) {
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
            callingForBothSamples(out, cfg, sample, region, v1, v2);
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
