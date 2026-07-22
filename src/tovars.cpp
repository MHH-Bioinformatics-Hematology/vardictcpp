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
                       : (double)v.highQualityReadsCount; // hi/lo signal-to-noise
            var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads))
                     + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads));
            var.duprate = vd.duprate();

            if (allele[0] == '+') {                    // insertion (rare in nonInsertion map)
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
            }
            var.vartype = classifyType(var.refallele, var.varallele);
            var.genotype = var.frequency < 0.5
                         ? var.refallele + "/" + var.varallele
                         : var.varallele + "/" + var.varallele;

            // Reference-context columns (flanks): 5 bp windows from the reference.
            for (int i = 5; i >= 1; --i) var.leftseq += ref.at(position - i);
            for (int i = 1; i <= 5; ++i) var.rightseq += ref.at(position + 1 + i - 1);

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
                for (int i = 5; i >= 1; --i) var.leftseq += ref.at(position - i);
                for (int i = 1; i <= 5; ++i) var.rightseq += ref.at(position + i);
                result.push_back(std::move(var));
            }
        }
    }
    return result;
}

} // namespace vardict
