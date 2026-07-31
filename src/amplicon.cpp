#include "amplicon.hpp"
#include "util.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <set>
#include <sstream>
#include <vector>

namespace vardict {

// ------------------------------------------------------------------------------------------------
// RegionBuilder.buildAmpRegions
// ------------------------------------------------------------------------------------------------
std::vector<std::vector<Region>> buildAmpRegions(const std::vector<std::string>& segRaws, const Config& cfg) {
    // AMP_BED_ROW_FORMAT = new BedRowFormat(chr 0, start 1, end 2, thickStart 6, thickEnd 7, gene 3)
    std::map<std::string, std::vector<Region>> tsegs; // chr -> regions (insertion order preserved below)
    std::vector<std::string> chrOrder;
    for (const std::string& line : segRaws) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> f;
        std::stringstream ss(line); std::string tok;
        while (std::getline(ss, tok, '\t')) f.push_back(tok);
        if (f.size() < 8) continue;
        Region r;
        r.chr = f[0];
        int start = std::stoi(f[1]);
        int end = std::stoi(f[2]);
        r.gene = f[3];
        int insertStart = std::stoi(f[6]);
        int insertEnd = std::stoi(f[7]);
        if (cfg.zeroBased && start < end) { start++; insertStart++; }
        r.start = start; r.end = end;
        r.insertStart = insertStart; r.insertEnd = insertEnd;
        if (tsegs.find(r.chr) == tsegs.end()) chrOrder.push_back(r.chr);
        tsegs[r.chr].push_back(r);
    }

    std::vector<std::vector<Region>> segs;
    segs.emplace_back();
    int previousEnd = -1;
    // Java iterates the HashMap entrySet (undefined order). Real amplicon BEDs are single-chromosome or
    // pre-grouped; we iterate chromosomes in first-seen order for determinism.
    for (const std::string& chr : chrOrder) {
        std::vector<Region>& chrRegions = tsegs[chr];
        std::stable_sort(chrRegions.begin(), chrRegions.end(),
                         [](const Region& a, const Region& b) { return a.insertStart < b.insertStart; });
        std::string previousChr;
        for (const Region& region : chrRegions) {
            if (previousEnd != -1 && (region.chr != previousChr || region.insertStart > previousEnd)) {
                segs.emplace_back();
            }
            segs.back().push_back(region);
            previousChr = region.chr;
            previousEnd = region.insertEnd;
        }
    }
    return segs;
}

// ------------------------------------------------------------------------------------------------
// Per-amplicon Vars (reuse the somatic candidate builder: same sorted, good-flagged, no-drop set)
// ------------------------------------------------------------------------------------------------
std::map<int, AmpVars> buildAmpVars(const Config& cfg, const Region& region,
                                    const VariationData& vd, Reference& ref) {
    std::map<int, AmpVars> result;
    std::vector<SomaticPosition> positions = callVariantsSomatic(cfg, region, vd, ref);
    for (SomaticPosition& sp : positions) {
        AmpVars av;
        av.variants = std::move(sp.variants);
        av.hasRef = true;
        auto covIt = vd.refCoverage.find(sp.position);
        av.refTotalCov = (covIt != vd.refCoverage.end()) ? covIt->second : 0;
        result[sp.position] = std::move(av);
    }
    return result;
}

// ------------------------------------------------------------------------------------------------
// Variant.adjComplex (trim shared prefix/suffix of a Complex variant's ref/alt)
// ------------------------------------------------------------------------------------------------
static void adjComplex(Variant& v) {
    std::string refAllele = v.refallele;
    std::string varAllele = v.varallele;
    if (!varAllele.empty() && varAllele[0] == '<') return;
    int n = 0;
    while ((int)refAllele.size() - n > 1 && (int)varAllele.size() - n > 1 &&
           refAllele[n] == varAllele[n]) n++;
    if (n > 0) {
        v.startPosition += n;
        v.refallele = substr(refAllele, n);
        v.varallele = substr(varAllele, n);
        v.leftseq += substr(refAllele, 0, n);
        v.leftseq = substr(v.leftseq, n);
    }
    refAllele = v.refallele;
    varAllele = v.varallele;
    n = 1;
    while ((int)refAllele.size() - n > 0 && (int)varAllele.size() - n > 0 &&
           substr(refAllele, -n, 1) == substr(varAllele, -n, 1)) n++;
    if (n > 1) {
        v.endPosition -= n - 1;
        v.refallele = substr(refAllele, 0, 1 - n);
        v.varallele = substr(varAllele, 0, 1 - n);
        v.rightseq = substr(refAllele, 1 - n, n - 1) + substr(v.rightseq, 0, 1 - n);
    }
}

// ------------------------------------------------------------------------------------------------
// AmpliconPostProcessModule helpers
// ------------------------------------------------------------------------------------------------
static int countVariantOnAmplicons(const Variant& vref, const std::map<int, std::vector<Variant>>& goodVariantsOnAmp) {
    int gvscnt = 0;
    for (const auto& [amp, variants] : goodVariantsOnAmp)
        for (const Variant& variant : variants)
            if (variant.refallele == vref.refallele && variant.varallele == vref.varallele) gvscnt++;
    return gvscnt;
}

static void fillVrefList(const std::vector<std::pair<Variant, std::string>>& gvs, std::vector<Variant>& vrefList) {
    for (const auto& gv : gvs) {
        bool added = false;
        for (const Variant& var : vrefList)
            if (var.varallele == gv.first.varallele && var.refallele == gv.first.refallele) added = true;
        if (!added) vrefList.push_back(gv.first);
    }
}

static bool isAmpBiasFlag(std::map<int, std::vector<Variant>>& goodVariantsOnAmp) {
    if (goodVariantsOnAmp.empty()) return false;
    std::vector<int> ampliconList;
    for (const auto& kv : goodVariantsOnAmp) ampliconList.push_back(kv.first);
    // std::map already sorts keys ascending (matches Collections.sort(ampliconList)).
    int ampliconLength = (int)ampliconList.size() - 1;
    auto tcovDesc = [](const Variant& a, const Variant& b) { return a.totalPosCoverage > b.totalPosCoverage; };
    for (int i = 0; i < ampliconLength; i++) {
        std::vector<Variant>& cur = goodVariantsOnAmp[ampliconList[i]];
        auto nit = goodVariantsOnAmp.find(ampliconList[i + 1]);
        if (nit == goodVariantsOnAmp.end() || cur.size() != nit->second.size()) return true;
        std::vector<Variant>& nxt = nit->second;
        std::stable_sort(cur.begin(), cur.end(), tcovDesc);
        std::stable_sort(nxt.begin(), nxt.end(), tcovDesc);
        for (size_t j = 0; j < cur.size(); j++)
            if (cur[j].descriptionString != nxt[j].descriptionString) return true;
    }
    return false;
}

// SimpleOutputVariant-style zero-or-fixed-precision formatter (== the amplicon 38-col format).
static std::string fmt(double v, const char* pat) {
    if (v == 0) return "0";
    char buf[64];
    std::snprintf(buf, sizeof(buf), pat, v);
    return buf;
}

// AmpliconOutputVariant.create_amplicon_variant_38columns
static void appendAmpliconRow(std::string& out, const Config& cfg, const Region& rg,
                              const Variant& variant, const std::string& seg,
                              int goodVariantsCount, int totalVariantsCount, int noCoverage, bool flag) {
    std::string af    = fmt(variant.frequency, "%.4f");
    std::string pmean = fmt(variant.pmean, "%.1f");
    std::string qmean = fmt(variant.qmean, "%.1f");
    std::string mq    = fmt(variant.mapq, "%.1f");
    std::string sn    = fmt(variant.qratio, "%.3f");
    std::string hiaf  = fmt(variant.hifreq, "%.4f");
    std::string exaf  = fmt(variant.extrafreq, "%.4f");
    std::string msi   = fmt(variant.msi, "%.3f");
    std::string nm    = variant.nm > 0 ? fmt(variant.nm, "%.1f") : std::string("0");
    const char* geno  = variant.genotype.empty() ? "0" : variant.genotype.c_str();
    const char* bias  = variant.bias.empty() ? "0" : variant.bias.c_str();

    // refallele/varallele/genotype/flanks of a large complex/insertion variant can be hundreds of bases;
    // a fixed buffer would truncate the line, dropping columns and the newline (merging the next row).
    std::vector<char> buf(512 + cfg.sample.size() + rg.gene.size() + rg.chr.size()
                          + variant.refallele.size() + variant.varallele.size() + variant.genotype.size()
                          + variant.leftseq.size() + variant.rightseq.size() + seg.size()
                          + variant.vartype.size() + variant.bias.size());
    std::snprintf(buf.data(), buf.size(),
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t"
        "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
        cfg.sample.c_str(),
        rg.gene.c_str(),
        rg.chr.c_str(),
        variant.startPosition,
        variant.endPosition,
        variant.refallele.c_str(),
        variant.varallele.c_str(),
        variant.totalPosCoverage,
        variant.varsCount,
        variant.refFwd, variant.refRev,
        variant.varFwd, variant.varRev,
        geno,
        af.c_str(),
        bias,
        pmean.c_str(),
        variant.pstd,
        qmean.c_str(),
        variant.qstd,
        mq.c_str(),
        sn.c_str(),
        hiaf.c_str(),
        exaf.c_str(),
        variant.shift3,
        msi.c_str(),
        variant.msint,
        nm.c_str(),
        variant.hicnt,
        variant.hicov,
        variant.leftseq.empty() ? "0" : variant.leftseq.c_str(),
        variant.rightseq.empty() ? "0" : variant.rightseq.c_str(),
        seg.c_str(),
        variant.vartype.c_str(),
        goodVariantsCount,
        totalVariantsCount,
        noCoverage,
        flag ? 1 : 0);
    out += buf.data();
}

// Reference-only row (variant == null): AmpliconOutputVariant(null, rg, null, null, position, 0, nocov, false).
static void appendAmpliconRefRow(std::string& out, const Config& cfg, const Region& rg,
                                 int position, int noCoverage) {
    std::string seg = rg.chr + ":" + std::to_string(position) + "-" + std::to_string(position);
    char buf[1024];
    // All variant-derived fields default to empty/zero; the printer emits "" for null string columns.
    std::snprintf(buf, sizeof(buf),
        "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t"
        "%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\n",
        cfg.sample.c_str(), rg.gene.c_str(), rg.chr.c_str(),
        position, position, "", "",
        0, 0, 0, 0, 0, 0,
        "", "0", "0;0", "0", 0, "0", 0,
        "0", "0", "0", "0", 0, "0", 0, "0", 0, 0,
        "0", "0", seg.c_str(), "", 0, 0, noCoverage, 0);
    out += buf;
}

// ------------------------------------------------------------------------------------------------
// AmpliconPostProcessModule.process for one segment
// ------------------------------------------------------------------------------------------------
void appendAmpliconSegment(std::string& out, const Config& cfg, const std::vector<Region>& regions,
                           const std::vector<std::map<int, AmpVars>>& vars) {
    if (regions.empty()) return;
    const Region& rg = regions.back(); // AmpliconMode passes currentRegion (last region of the segment)

    // ampliconsOnPositions: position -> list of (ampliconNumber). ampliconNumber indexes regions/vars.
    std::map<int, std::vector<int>> pos;
    for (int a = 0; a < (int)regions.size(); ++a)
        for (int p = regions[a].insertStart; p <= regions[a].insertEnd; ++p)
            pos[p].push_back(a);

    for (const auto& [position, ampliconRegions] : pos) {
        std::vector<std::pair<Variant, std::string>> gvs; // good variants (variant, "chr:start-end")
        std::vector<Variant> ref;                          // reference variants (only coverage matters)
        std::vector<Variant> vrefList;
        std::set<std::string> goodmap;
        std::vector<int> vcovs;
        std::map<int, std::vector<Variant>> goodVariantsOnAmp;
        int nocov = 0, maxcov = 0; double maxaf = 0;

        for (int ampliconNumber : ampliconRegions) {
            const Region& reg = regions[ampliconNumber];
            const std::string& chr = reg.chr;
            int start = reg.start, end = reg.end;

            auto vit = vars[ampliconNumber].find(position);
            const AmpVars* vtmp = (vit == vars[ampliconNumber].end()) ? nullptr : &vit->second;
            const std::vector<Variant>* variantsOnAmplicon = vtmp ? &vtmp->variants : nullptr;

            if (variantsOnAmplicon && !variantsOnAmplicon->empty()) {
                std::vector<Variant> goodVars;
                for (const Variant& tv : *variantsOnAmplicon) {
                    vcovs.push_back(tv.totalPosCoverage);
                    if (tv.totalPosCoverage > maxcov) maxcov = tv.totalPosCoverage;
                    if (tv.good) {
                        gvs.push_back({tv, chr + ":" + std::to_string(start) + "-" + std::to_string(end)});
                        goodVars.push_back(tv);
                        goodVariantsOnAmp[ampliconNumber] = goodVars;
                        if (tv.frequency > maxaf) maxaf = tv.frequency;
                        goodmap.insert(std::to_string(ampliconNumber) + "-" + tv.refallele + "-" + tv.varallele);
                    }
                }
            } else if (vtmp && vtmp->hasRef) {
                vcovs.push_back(vtmp->refTotalCov);
            } else {
                vcovs.push_back(0);
            }
            if (vtmp && vtmp->hasRef) { Variant rv; rv.totalPosCoverage = vtmp->refTotalCov; ref.push_back(rv); }
        }

        for (int t : vcovs) if (t < maxcov / (double)50) nocov++;

        if (gvs.size() > 1)
            std::stable_sort(gvs.begin(), gvs.end(),
                             [](const auto& a, const auto& b) { return b.first.frequency < a.first.frequency; });
        if (ref.size() > 1)
            std::stable_sort(ref.begin(), ref.end(),
                             [](const Variant& a, const Variant& b) { return a.totalPosCoverage > b.totalPosCoverage; });

        if (gvs.empty()) { // only reference
            if (cfg.doPileup) {
                if (!ref.empty()) {
                    vrefList.push_back(ref[0]);
                } else {
                    appendAmpliconRefRow(out, cfg, rg, position, nocov);
                    continue;
                }
            } else {
                continue;
            }
        } else {
            fillVrefList(gvs, vrefList);
        }
        bool flag = isAmpBiasFlag(goodVariantsOnAmp);

        for (size_t i = 0; i < vrefList.size(); ++i) {
            Variant vref = vrefList[i];
            std::vector<std::pair<Variant, std::string>> goodVariants = gvs;
            if (flag) { // different good variants across amplicons: re-collect this variant's amplicons
                const std::string& gdnt = gvs[0].first.descriptionString;
                std::vector<std::pair<Variant, std::string>> gcnt;
                for (int amp : ampliconRegions) {
                    auto vit = vars[amp].find(position);
                    const AmpVars* vtmp = (vit == vars[amp].end()) ? nullptr : &vit->second;
                    if (!vtmp) continue;
                    for (const Variant& variant : vtmp->variants) {
                        if (variant.descriptionString == gdnt && variant.good) {
                            const Region& r2 = regions[amp];
                            gcnt.push_back({variant, r2.chr + ":" + std::to_string(r2.start) + "-" + std::to_string(r2.end)});
                            break;
                        }
                    }
                }
                if (gcnt.size() == gvs.size()) flag = false;
                std::stable_sort(gcnt.begin(), gcnt.end(),
                                 [](const auto& a, const auto& b) { return b.first.frequency < a.first.frequency; });
                goodVariants = gcnt;
            }
            int initialGvscnt = countVariantOnAmplicons(vref, goodVariantsOnAmp);
            int currentGvscnt = initialGvscnt;
            std::vector<std::pair<Variant, std::string>> badVariants;
            if (initialGvscnt != (int)ampliconRegions.size() || flag) {
                for (int amp : ampliconRegions) {
                    const Region& reg = regions[amp];
                    if (goodmap.count(std::to_string(amp) + "-" + vref.refallele + "-" + vref.varallele)) continue;
                    if (cfg.doPileup && vref.refallele == vref.varallele) continue;
                    if (vref.startPosition >= reg.insertStart && vref.endPosition <= reg.insertEnd) {
                        std::string regStr = reg.chr + ":" + std::to_string(reg.start) + "-" + std::to_string(reg.end);
                        auto vit = vars[amp].find(position);
                        const AmpVars* vtmp = (vit == vars[amp].end()) ? nullptr : &vit->second;
                        if (vtmp && !vtmp->variants.empty()) {
                            badVariants.push_back({vtmp->variants[0], regStr});
                        } else if (vtmp && vtmp->hasRef) {
                            Variant rv; rv.totalPosCoverage = vtmp->refTotalCov;
                            badVariants.push_back({rv, regStr});
                        } else {
                            badVariants.push_back({Variant(), regStr});
                        }
                    } else if ((vref.startPosition < reg.insertEnd && reg.insertEnd < vref.endPosition) ||
                               (vref.startPosition < reg.insertStart && reg.insertStart < vref.endPosition)) {
                        if (currentGvscnt > 1) currentGvscnt--;
                    }
                }
            }
            if (flag && currentGvscnt < initialGvscnt) flag = false;
            if (vref.vartype == "Complex") adjComplex(vref);
            std::string seg = !goodVariants.empty()
                              ? goodVariants[0].second
                              : rg.chr + ":" + std::to_string(position) + "-" + std::to_string(position);
            int totalVariantsCount = currentGvscnt + (int)badVariants.size();
            appendAmpliconRow(out, cfg, rg, vref, seg, currentGvscnt, totalVariantsCount, nocov, flag);
        }
    }
}

void printAmpliconHeader(std::FILE* out) {
    std::fprintf(out,
        "Sample\tGene\tChr\tStart\tEnd\tRef\tAlt\tDepth\tAltDepth\tRefFwdReads\tRefRevReads\t"
        "AltFwdReads\tAltRevReads\tGenotype\tAF\tBias\tPMean\tPStd\tQMean\tQStd\tMQ\tSig_Noise\t"
        "HiAF\tExtraAF\tshift3\tMSI\tMSI_NT\tNM\tHiCnt\tHiCov\t5pFlankSeq\t3pFlankSeq\tSeg\t"
        "VarType\tGoodVarCount\tTotalVarCount\tNocov\tAmpflag\n");
}

} // namespace vardict
