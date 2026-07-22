#include "tovars.hpp"
#include <algorithm>
#include <cmath>

namespace vardict {

// VarDict strand-bias flag (data/VariationUtils.strandBias). 0 = one strand only / too few,
// 1 = biased, 2 = both strands well represented.
static int strandBias(int fwd, int rev, int minBiasReads) {
    if (fwd + rev <= 12) {
        return (fwd > 0 && rev > 0) ? 2 : 0;
    }
    double tot = fwd + rev;
    bool ok = (fwd / tot >= 0.01) && (rev / tot >= 0.01) && fwd >= minBiasReads && rev >= minBiasReads;
    return ok ? 2 : 1;
}

// Port of variations/Variant.java isGoodVar for the simple single-sample path. Reference-allele
// stats (hicnt, mean mapping quality) are passed in from the position's ref accumulator. MSI columns
// default to 0 until findMSI is ported, so the two MSI gates are inactive (matches a no-MSI position).
static bool isGoodVar(const Config& c, const Variant& v, int refHicnt, double refMeanMapq) {
    if (v.refallele.empty()) return false;
    if (v.frequency < c.freq || v.hicnt < c.minReads ||
        v.pmean < c.readPosFilter || v.qmean < c.goodq) {
        return false;
    }
    if (refHicnt > c.minReads && v.frequency < 0.25) {
        double d = v.mapq + (double)v.refallele.size() + (double)v.varallele.size();
        double f = (1 + d) / (refMeanMapq + 1);
        if ((d - 2 < 5 && refMeanMapq > 20) || f < 0.25) return false;
    }
    // (Deletion/splice gate omitted: splice set not tracked yet.)
    if (v.qratio < c.qratio) return false;
    if (v.frequency > 0.30) return true;
    if (v.mapq < c.mapqMin) return false;
    if (v.msi >= 15 && v.frequency <= c.monomerMsiFrequency && v.msint == 1) return false;
    if (v.msi >= 12 && v.frequency <= c.nonMonomerMsiFrequency && v.msint > 1) return false;
    if (v.bias == "2;1" && v.frequency < 0.20) {
        if (v.vartype == "SNV" || (v.refallele.size() < 3 && v.varallele.size() < 3)) return false;
    }
    return true;
}

// Port of ToVarsBuilder.findMSI. tseq1 = left ref window ending at/after the variant, tseq2 = right
// ref window. Returns {msi count, shift3, microsatellite-unit length}. Repeat counting is done
// directly (alleles are ACGTN) rather than via regex, matching ((unit)+)$ on tseq1 and ^((unit)+) on tseq2.
struct MSIResult { double msi; int shift3; int msintLen; };
static MSIResult findMSI(const std::string& tseq1, const std::string& tseq2) {
    int nmsi = 1;
    double msicnt = 0;
    std::string maxmsi;
    while (nmsi <= (int)tseq1.size() && nmsi <= 6) {
        std::string unit = tseq1.substr(tseq1.size() - nmsi); // substr(tseq1, -nmsi)
        // trailing repeats of `unit` in tseq1
        int t1 = 0;
        while ((int)tseq1.size() - (t1 + 1) * nmsi >= 0 &&
               tseq1.compare(tseq1.size() - (t1 + 1) * nmsi, nmsi, unit) == 0) t1++;
        double curmsi = (double)(t1 * nmsi) / nmsi;
        // leading repeats of `unit` in tseq2
        int t2 = 0;
        while ((t2 + 1) * nmsi <= (int)tseq2.size() &&
               tseq2.compare(t2 * nmsi, nmsi, unit) == 0) t2++;
        curmsi += (double)(t2 * nmsi) / nmsi;
        if (curmsi > msicnt) { maxmsi = unit; msicnt = curmsi; }
        nmsi++;
    }
    std::string tseq = tseq1 + tseq2;
    int shift3 = 0;
    while (shift3 < (int)tseq2.size() && tseq[shift3] == tseq2[shift3]) shift3++;
    return { msicnt, shift3, (int)maxmsi.size() };
}

static std::string classifyType(const std::string& ref, const std::string& alt) {
    if (alt.size() == 1 && ref.size() == 1) return "SNV";
    if (!alt.empty() && alt[0] == '+') return "Insertion";
    if (!alt.empty() && alt[0] == '-') return "Deletion";
    if (ref.size() == alt.size()) return alt.size() == 1 ? "SNV" : "Complex";
    return "Complex";
}

std::vector<Variant> callVariants(const Config& cfg, const Region& region,
                                  const VariationData& vd, Reference& ref) {
    std::vector<Variant> result;
    const double freq = cfg.freq;

    for (const auto& [position, alleleMap] : vd.nonInsertionVariants) {
        if (position < region.start || position > region.end) continue;
        auto covIt = vd.refCoverage.find(position);
        if (covIt == vd.refCoverage.end() || covIt->second == 0) continue;
        int totalCov = covIt->second;
        char refBase = ref.at(position);

        // Reference-allele accumulator (allele == ref base).
        const Variation* refVar = nullptr;
        auto rit = alleleMap.find(std::string(1, refBase));
        if (rit != alleleMap.end()) refVar = &rit->second;
        int refHicnt = refVar ? refVar->highQualityReadsCount : 0;
        double refMeanMapq = (refVar && refVar->varsCount) ? refVar->meanMappingQuality / refVar->varsCount : 0;

        // High-quality coverage = sum of hi-qual reads across alleles at this position.
        int hicov = 0;
        for (const auto& [al, v] : alleleMap) hicov += v.highQualityReadsCount;

        for (const auto& [allele, v] : alleleMap) {
            if (allele.size() == 1 && allele[0] == refBase) continue; // skip pure reference
            if (v.varsCount < cfg.minReads) continue;
            double af = totalCov > 0 ? (double)v.varsCount / (double)totalCov : 0.0;
            if (!cfg.doPileup && af <= freq) continue; // -f filter (freq=0 keeps af>0)

            Variant var;
            var.startPosition = position;
            var.endPosition = position;
            var.totalPosCoverage = totalCov;
            var.varsCount = v.varsCount;
            var.varFwd = v.varsCountOnForward;
            var.varRev = v.varsCountOnReverse;
            if (refVar) { var.refFwd = refVar->varsCountOnForward; var.refRev = refVar->varsCountOnReverse; }
            var.frequency = af;
            var.pmean = v.varsCount ? v.meanPosition / v.varsCount : 0;
            var.qmean = v.varsCount ? v.meanQuality / v.varsCount : 0;
            var.mapq  = v.varsCount ? v.meanMappingQuality / v.varsCount : 0;
            var.nm    = v.varsCount ? v.numberOfMismatches / v.varsCount : 0;
            var.pstd  = v.pstd ? 1 : 0;
            var.qstd  = v.qstd ? 1 : 0;
            var.hicnt = v.highQualityReadsCount;
            var.hicov = hicov;
            var.hifreq = hicov > 0 ? (double)v.highQualityReadsCount / hicov : 0;
            var.qratio = v.lowQualityReadsCount > 0
                       ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                       : (double)v.highQualityReadsCount / 0.5; // hi/lo signal-to-noise
            var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads))
                     + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads));
            var.duprate = vd.duprate();

            if (allele.find('&') != std::string::npos) {   // MNV / complex (e.g. "A&CG" -> ACG)
                std::string va;
                for (char ch : allele) if (ch != '&') va += ch;
                var.varallele = va;
                std::string ra;
                for (int k = 0; k < (int)va.size(); ++k) ra += ref.at(position + k);
                var.refallele = ra;
                var.endPosition = position + (int)va.size() - 1;
                // MSI for MNP/complex uses the same SNV/MNP reference window as findMSI.
                std::string t1, t2;
                for (int q = position - 30; q <= position + 1; ++q) if (q >= 1) t1 += ref.at(q);
                for (int q = position + 2; q <= position + 70; ++q) t2 += ref.at(q);
                MSIResult m = findMSI(t1, t2);
                var.msi = m.msi; var.shift3 = m.shift3; var.msint = m.msintLen;
            } else if (allele[0] == '+') {              // insertion (rare in nonInsertion map)
                var.refallele = std::string(1, refBase);
                var.varallele = std::string(1, refBase) + allele.substr(1);
            } else if (allele[0] == '-') {             // deletion signature "-N"
                int dl = std::stoi(allele.substr(1));
                std::string delseq;
                for (int i = 0; i < dl; ++i) delseq += ref.at(position + i);
                var.refallele = std::string(1, refBase) + delseq;
                var.varallele = std::string(1, refBase);
                var.endPosition = position + dl;
            } else {                                    // SNV
                var.refallele = std::string(1, refBase);
                var.varallele = allele;
                // MSI adjustment for SNV/MNP (ToVarsBuilder: findMSI on ref[p-30..p+1], ref[p+2..p+70]).
                std::string tseq1, tseq2;
                for (int q = position - 30; q <= position + 1; ++q) if (q >= 1) tseq1 += ref.at(q);
                for (int q = position + 2; q <= position + 70; ++q) tseq2 += ref.at(q);
                MSIResult m = findMSI(tseq1, tseq2);
                var.msi = m.msi; var.shift3 = m.shift3; var.msint = m.msintLen;
            }
            var.vartype = classifyType(var.refallele, var.varallele);
            var.genotype = var.frequency < 0.5
                         ? var.refallele + "/" + var.varallele
                         : var.varallele + "/" + var.varallele;

            // Reference-context flanks: 20 bp windows (ToVarsBuilder REF_20_BASES).
            for (int i = 20; i >= 1; --i) if (position - i >= 1) var.leftseq += ref.at(position - i);
            for (int i = 1; i <= 20; ++i) var.rightseq += ref.at(var.endPosition + i);

            if (!cfg.doPileup && !isGoodVar(cfg, var, refHicnt, refMeanMapq)) continue;
            result.push_back(std::move(var));
        }

        // Insertions anchored at this position.
        auto insIt = vd.insertionVariants.find(position);
        if (insIt != vd.insertionVariants.end()) {
            for (const auto& [allele, v] : insIt->second) {
                if (v.varsCount < cfg.minReads) continue;
                double af = totalCov > 0 ? (double)v.varsCount / (double)totalCov : 0.0;
                if (!cfg.doPileup && af <= freq) continue;
                Variant var;
                var.startPosition = position; var.endPosition = position;
                var.totalPosCoverage = totalCov; var.varsCount = v.varsCount;
                var.varFwd = v.varsCountOnForward; var.varRev = v.varsCountOnReverse;
                if (refVar) { var.refFwd = refVar->varsCountOnForward; var.refRev = refVar->varsCountOnReverse; }
                var.frequency = af;
                var.pmean = v.varsCount ? v.meanPosition / v.varsCount : 0;
                var.qmean = v.varsCount ? v.meanQuality / v.varsCount : 0;
                var.mapq  = v.varsCount ? v.meanMappingQuality / v.varsCount : 0;
                var.nm    = v.varsCount ? v.numberOfMismatches / v.varsCount : 0;
                var.hicnt = v.highQualityReadsCount; var.hicov = hicov;
                var.refallele = std::string(1, refBase);
                var.varallele = std::string(1, refBase) + allele.substr(1);
                var.vartype = "Insertion";
                var.genotype = var.refallele + "/" + var.varallele;
                var.duprate = vd.duprate();
                var.qratio = v.lowQualityReadsCount > 0
                           ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                           : (double)v.highQualityReadsCount / 0.5;
                var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads))
                         + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads));
                for (int i = 20; i >= 1; --i) if (position - i >= 1) var.leftseq += ref.at(position - i);
                for (int i = 1; i <= 20; ++i) var.rightseq += ref.at(position + i);
                if (!cfg.doPileup && !isGoodVar(cfg, var, refHicnt, refMeanMapq)) continue;
                result.push_back(std::move(var));
            }
        }
    }
    return result;
}

} // namespace vardict
