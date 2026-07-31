#include "tovars.hpp"
#include "util.hpp"
#include <algorithm>
#include <cmath>

namespace vardict {

// VarDict strand-bias flag (data/VariationUtils.strandBias). 0 = one strand only / too few,
// 1 = biased, 2 = both strands well represented.
static int strandBias(int fwd, int rev, int minBiasReads, double bias) {
    if (fwd + rev <= 12) {
        return (fwd > 0 && rev > 0) ? 2 : 0;
    }
    double tot = fwd + rev;
    bool ok = (fwd / tot >= bias) && (rev / tot >= bias) && fwd >= minBiasReads && rev >= minBiasReads;
    return ok ? 2 : 1;
}

// Ports ToVarsBuilder's complex-descriptor decoding (lines 765-846): AMP_ATGC "&([ATGC]+)",
// HASH_GROUP_CARET_GROUP "#(.+)\\^(.+)" and CARET_ATGNC "\\^([ATGNC]+)". Turns a raw CigarParser
// complex description (e.g. insertion "+A&CT", deletion "-2&AC^G#TT") into the realized
// refallele/varallele plus adjusted start/end positions, mirroring VarDictJava byte-for-byte.
// `descriptionString` is the original allele; refallele/varallele are seeded by the caller.
static void applyComplexGrammar(const std::string& descriptionString, Reference& ref,
                                std::string& refallele, std::string& varallele,
                                int& startPos, int& endPos, std::string& genotype1current) {
    auto joinRefLocal = [&](int from, int to) {
        std::string r; for (int i = from; i <= to; ++i) if (ref.has(i)) r += ref.at(i); return r;
    };
    // First "&([ATGC]+)" group in s (>=1 ATGC after '&'); empty string == no match.
    auto ampGroup = [](const std::string& s) -> std::string {
        auto amp = s.find('&'); if (amp == std::string::npos) return "";
        std::string g;
        for (size_t i = amp + 1; i < s.size(); ++i) { char c = s[i]; if (c=='A'||c=='T'||c=='G'||c=='C') g += c; else break; }
        return g;
    };
    auto eraseFirst = [](std::string& s, char c) { auto p = s.find(c); if (p != std::string::npos) s.erase(p, 1); };

    // AMP_ATGC: variant followed by matched reference sequence.
    std::string extra = ampGroup(descriptionString);
    if (!extra.empty()) {
        eraseFirst(varallele, '&');
        std::string tch = joinRefLocal(endPos + 1, endPos + (int)extra.size());
        refallele += tch; genotype1current += tch; endPos += (int)extra.size();
        std::string vextra = ampGroup(varallele);
        if (!vextra.empty()) {
            eraseFirst(varallele, '&');
            std::string tch2 = joinRefLocal(endPos + 1, endPos + (int)vextra.size());
            refallele += tch2; genotype1current += tch2; endPos += (int)vextra.size();
        }
        if (!descriptionString.empty() && descriptionString[0] == '+') {
            if (!refallele.empty()) refallele = refallele.substr(1);
            if (!varallele.empty()) varallele = varallele.substr(1);
            startPos++;
        }
    }

    // HASH_GROUP_CARET_GROUP "#(.+)\\^(.+)": short matched sequence + indel tail. Greedy group(1)
    // runs from the first '#' to the LAST '^'.
    {
        auto hp = descriptionString.find('#');
        auto cp = descriptionString.rfind('^');
        if (hp != std::string::npos && cp != std::string::npos && cp > hp + 1 && cp + 1 < descriptionString.size()) {
            std::string matchedSequence = descriptionString.substr(hp + 1, cp - (hp + 1));
            std::string tail = descriptionString.substr(cp + 1);
            endPos += (int)matchedSequence.size();
            refallele += joinRefLocal(endPos - (int)matchedSequence.size() + 1, endPos);
            int d = 0; size_t z = 0; while (z < tail.size() && isdigit((unsigned char)tail[z])) { d = d*10 + (tail[z]-'0'); z++; }
            if (z > 0) { refallele += joinRefLocal(endPos + 1, endPos + d); endPos += d; }
            eraseFirst(varallele, '#');
            auto vc = varallele.find('^');
            if (vc != std::string::npos) { size_t e = vc + 1; while (e < varallele.size() && isdigit((unsigned char)varallele[e])) e++; varallele.erase(vc, e - vc); }
        }
    }

    // CARET_ATGNC "\\^([ATGNC]+)": deletion followed directly by insertion.
    {
        auto cp = descriptionString.find('^');
        if (cp != std::string::npos && cp + 1 < descriptionString.size()) {
            char c = descriptionString[cp + 1];
            if (c=='A'||c=='T'||c=='G'||c=='N'||c=='C') { auto vc = varallele.find('^'); if (vc != std::string::npos) varallele.erase(vc, 1); }
        }
    }
}

// Variant.adjComplex (SimplePostProcessModule): trim the shared 5' prefix and 3' suffix of a Complex
// variant's ref/alt (keeping >= 1 base each side), shifting start/end and the flanking sequences.
// Applied to Complex variants only, exactly as VarDictJava's post-processing does.
static void adjComplexVar(Variant& v) {
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

// Port of variations/Variant.java isGoodVar for the simple single-sample path. Reference-allele
// stats (hicnt, mean mapping quality) are passed in from the position's ref accumulator. MSI columns
// default to 0 until findMSI is ported, so the two MSI gates are inactive (matches a no-MSI position).
static bool isGoodVar(const Config& c, const Variant& v, int refHicnt, double refMeanMapq,
                      const std::set<std::string>& splice) {
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
    // A deletion whose coordinates exactly match a recorded splice junction is an intron, not a variant.
    if (v.vartype == "Deletion" &&
        splice.count(std::to_string(v.startPosition) + "-" + std::to_string(v.endPosition))) return false;
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
static MSIResult findMSI(const std::string& tseq1, const std::string& tseq2, const std::string& left = "") {
    int nmsi = 1;
    double msicnt = 0;
    std::string maxmsi;
    while (nmsi <= (int)tseq1.size() && nmsi <= 6) {
        std::string unit = tseq1.substr(tseq1.size() - nmsi); // substr(tseq1, -nmsi)
        // trailing repeats of `unit`; when `left` is given, count over (left + tseq1) as VarDict does.
        std::string b = left.empty() ? tseq1 : (left + tseq1);
        int t1 = 0;
        while ((int)b.size() - (t1 + 1) * nmsi >= 0 &&
               b.compare(b.size() - (t1 + 1) * nmsi, nmsi, unit) == 0) t1++;
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

// Round to 4 decimals with round-half-to-even (matches Java Utils.roundHalfEven("0.0000", x));
// %.4f uses the default IEEE round-to-nearest-even. VarDict stores the *rounded* frequency in
// createVariant BEFORE the position-level `maxfreq <= freq` filter, so a variant at 3/299 = 0.010033
// rounds to 0.0100 and is dropped; comparing the unrounded value would wrongly keep it.
static double roundN(double x, int n) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.*f", n, x);
    return std::atof(buf);
}
static double round4(double x) { return roundN(x, 4); }

// isGoodVar in Java tests the *stored* (rounded) pmean/qmean/mapq/msi and referenceVar.meanMappingQuality
// (ToVarsBuilder rounds HALF_EVEN before storing: pmean/qmean/mapq to 1 dp, msi to 3 dp). Round the
// filter probe the same way so boundary rows (e.g. qmean 22.46 -> 22.5) match Java's keep/drop decision.
static void roundProbeForFilter(Variant& probe) {
    probe.frequency = round4(probe.frequency);
    probe.pmean = roundN(probe.pmean, 1);
    probe.qmean = roundN(probe.qmean, 1);
    probe.mapq  = roundN(probe.mapq, 1);
    probe.msi   = roundN(probe.msi, 3);
}

// Port of Variant.varType() (classify by realized ref/alt alleles, not the raw description).
static std::string classifyType(const std::string& ref, const std::string& alt) {
    if (ref == alt && ref.size() == 1) return "";
    if (ref.size() == 1 && alt.size() == 1) return "SNV";
    if (!alt.empty() && (alt.front() == '<' || alt[0] == '+' || alt[0] == '-')) {
        if (alt[0] == '+') return "Insertion";
        if (alt[0] == '-') return "Deletion";
        // <DEL>/<DUP>/<INV>: Variant.varType extracts the 3-letter type between the angle brackets.
        if (alt.size() >= 3 && alt.back() == '>') return alt.substr(1, alt.size() - 2);
        return alt;
    }
    if (ref.empty() || alt.empty()) return "Complex";
    if (ref[0] != alt[0]) return "Complex";
    if (ref.size() == 1 && alt.size() > 1 && alt.compare(0, ref.size(), ref) == 0) return "Insertion";
    if (ref.size() > 1 && alt.size() == 1 && ref.compare(0, alt.size(), alt) == 0) return "Deletion";
    return "Complex";
}

std::vector<Variant> callVariants(const Config& cfg, const Region& region,
                                  const VariationData& vd, Reference& ref) {
    std::vector<Variant> result;
    const double freq = cfg.freq;

    // nonInsertionVariants is an unordered_map (fast per-base insertion); emit in ascending position.
    std::vector<int> sortedPositions;
    sortedPositions.reserve(vd.nonInsertionVariants.size());
    for (const auto& kv : vd.nonInsertionVariants) sortedPositions.push_back(kv.first);
    std::sort(sortedPositions.begin(), sortedPositions.end());
    for (int position : sortedPositions) {
        const VarMap& alleleMap = vd.nonInsertionVariants.at(position);
        // A position carrying a structural-variant marker is emitted even outside the region window
        // (ToVarsBuilder skips the region-bounds check when varsAtCurPosition.sv != null).
        auto svIt = vd.svInfoAt.find(position);
        bool isSVpos = svIt != vd.svInfoAt.end();
        if (!isSVpos && (position < region.start || position > region.end)) continue;
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

        // genotype1 = dominant description string, shared by every variant row at this position
        // (ToVarsBuilder.collectReferenceVariants). If the reference variant has frequency >= -f it
        // is the reference base; otherwise it is the highest-ranked non-reference variant, where the
        // rank key is roundHalfEven("0.0", meanQuality) * varsCount descending, tie broken by the
        // description string ascending (sortVariants). Candidates include non-insertion variants AND
        // insertions (createVariant + createInsertion feed the same sorted list).
        std::string positionGenotype1;
        {
            double refFreqG = (refVar && totalCov > 0) ? (double)refVar->varsCount / totalCov : 0;
            if (refVar && refFreqG >= freq) {
                positionGenotype1 = std::string(1, refBase);
            } else {
                const std::string* best = nullptr;
                double bestKey = 0;
                auto consider = [&](const std::string& al, const Variation& vv) {
                    if (al.size() == 1 && al[0] == refBase) return; // reference variant, not in list
                    if (vv.varsCount == 0) return;
                    double q = roundHalfEven("0.0", vv.meanQuality / (double)vv.varsCount);
                    double key = q * (double)vv.varsCount;
                    if (best == nullptr || key > bestKey || (key == bestKey && al < *best)) {
                        best = &al; bestKey = key;
                    }
                };
                for (const auto& [al, vv] : alleleMap) consider(al, vv);
                auto insG = vd.insertionVariants.find(position);
                if (insG != vd.insertionVariants.end())
                    for (const auto& [al, vv] : insG->second) consider(al, vv);
                if (best) positionGenotype1 = *best;
                else positionGenotype1 = std::string(1, refBase);
            }
            // '+' handling (collectReferenceVariants): plain insertion -> "+<insertedLength>".
            if (!positionGenotype1.empty() && positionGenotype1[0] == '+' &&
                positionGenotype1.find("<dup") == std::string::npos) {
                positionGenotype1 = "+" + std::to_string((int)positionGenotype1.size() - 1);
            }
        }

        // Per-position frequency gate (ToVarsBuilder collectVarsAtPosition + the `maxfreq <= -f`
        // drop): the WHOLE position is kept iff the maximum rounded AF over all non-reference
        // variants (non-insertions AND insertions) exceeds -f. A sub-threshold variant (e.g. a SNV at
        // exactly 0.0100) therefore survives when a sibling variant at the same position clears -f -
        // this is NOT a per-variant filter. maxfreq mirrors createVariant/createInsertion's ttcov.
        {
            double maxfreq = 0.0;
            for (const auto& [al, vv] : alleleMap) {
                if (vv.varsCount == 0) continue;
                if (al.size() == 1 && al[0] == refBase) continue;   // reference variant excluded
                int ttcov = totalCov;
                if (vv.varsCount > totalCov && vv.extracnt != 0 && vv.varsCount - totalCov < vv.extracnt) ttcov = vv.varsCount;
                double f = ttcov > 0 ? round4((double)vv.varsCount / ttcov) : 0.0;
                if (f > maxfreq) maxfreq = f;
            }
            auto insM = vd.insertionVariants.find(position);
            if (insM != vd.insertionVariants.end()) {
                int runningCov = totalCov;                            // createInsertion mutates totalPosCoverage across insertions
                for (const auto& [al, vv] : insM->second) {           // std::map iterates in sorted key order
                    if (vv.varsCount == 0) continue;
                    if (al.find('&') != std::string::npos) {          // '&' insertion re-bases coverage on position+1
                        auto c1 = vd.refCoverage.find(position + 1);
                        if (c1 != vd.refCoverage.end()) runningCov = c1->second;
                    }
                    int ttcov = runningCov;
                    if (vv.varsCount > runningCov && vv.extracnt != 0 && vv.varsCount - runningCov < vv.extracnt) ttcov = vv.varsCount;
                    if (ttcov < vv.varsCount) {
                        ttcov = vv.varsCount;
                        auto c1 = vd.refCoverage.find(position + 1);
                        if (c1 != vd.refCoverage.end() && ttcov < c1->second - vv.varsCount) ttcov = c1->second;
                        runningCov = ttcov;
                    }
                    double f = ttcov > 0 ? round4((double)vv.varsCount / ttcov) : 0.0;
                    if (f > maxfreq) maxfreq = f;
                }
            }
            if (!cfg.doPileup && freq > 0 && maxfreq <= freq) continue;   // drop the whole position
        }

        // createInsertion coverage reconciliation (Findings 2a,2c,2d) + collectReferenceVariants pos+1
        // swap (2b), computed read-only (vd is const): running over insertions in sorted order, the
        // position coverage is bumped for '&'/dominant insertions, and the dominant insertion's fwd/rev
        // are subtracted from the pos+1 reference variant. The final coverage is the Depth of every
        // variant here, and when it exceeds refCoverage[position] the RefFwd/RefRev are re-sourced from
        // that (reduced) pos+1 reference variant.
        int finalTotalCov = totalCov;
        int refFwdOut = refVar ? refVar->varsCountOnForward : 0;
        int refRevOut = refVar ? refVar->varsCountOnReverse : 0;
        // Per-insertion ttcov (ToVarsBuilder.createInsertion): the frequency/extraFrequency denominator
        // for each insertion, which can EXCEED the position Depth (the running totalPosCoverage is only
        // bumped in the `ttcov < varsCount` branch, but the extracnt branch raises ttcov without it).
        std::map<std::string,int> insTtcov;
        {
            auto c1it = vd.refCoverage.find(position + 1);
            int cov1 = (c1it != vd.refCoverage.end()) ? c1it->second : -1;
            int subFwd = 0, subRev = 0;
            auto insR = vd.insertionVariants.find(position);
            if (insR != vd.insertionVariants.end()) {
                int runningCov = totalCov;
                for (const auto& [al, vv] : insR->second) {     // std::map sorted == Java Collections.sort
                    if (al.find('&') != std::string::npos && cov1 >= 0) runningCov = cov1;
                    int ttcov = runningCov;
                    if (vv.varsCount > runningCov && vv.extracnt != 0 && vv.varsCount - runningCov < vv.extracnt) ttcov = vv.varsCount;
                    if (ttcov < vv.varsCount) {
                        ttcov = vv.varsCount;
                        if (cov1 >= 0 && ttcov < cov1 - vv.varsCount) {
                            ttcov = cov1;
                            subFwd += vv.varsCountOnForward; subRev += vv.varsCountOnReverse;
                        }
                        runningCov = ttcov;
                    }
                    insTtcov[al] = ttcov;
                }
                finalTotalCov = runningCov;
            }
            if (finalTotalCov > totalCov) {   // totalCov == refCoverage[position]
                auto p1 = vd.nonInsertionVariants.find(position + 1);
                if (p1 != vd.nonInsertionVariants.end() && ref.has(position + 1)) {
                    auto rv1 = p1->second.find(std::string(1, ref.at(position + 1)));
                    if (rv1 != p1->second.end()) {
                        refFwdOut = rv1->second.varsCountOnForward - subFwd;
                        refRevOut = rv1->second.varsCountOnReverse - subRev;
                    }
                }
            }
        }

        for (const auto& [allele, v] : alleleMap) {
            if (allele.size() == 1 && allele[0] == refBase) continue; // skip pure reference
            if (v.varsCount < cfg.minReads) continue;
            double af = totalCov > 0 ? (double)v.varsCount / (double)totalCov : 0.0;

            Variant var;
            var.startPosition = position;
            var.endPosition = position;
            var.totalPosCoverage = finalTotalCov;
            var.varsCount = v.varsCount;
            var.varFwd = v.varsCountOnForward;
            var.varRev = v.varsCountOnReverse;
            var.refFwd = refFwdOut; var.refRev = refRevOut;
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
            var.extrafreq = (v.extracnt != 0 && totalCov > 0) ? (double)v.extracnt / totalCov : 0;
            var.qratio = v.lowQualityReadsCount > 0
                       ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                       : (double)v.highQualityReadsCount / 0.5; // hi/lo signal-to-noise
            var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias))
                     + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias));
            var.duprate = vd.duprate();

            if (allele.find('&') != std::string::npos && allele[0] != '-' && allele[0] != '+') {   // MNV / complex (e.g. "A&CG" -> ACG)
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
            } else if (allele[0] == '-' && allele.find("<inv") != std::string::npos) {
                // Split-read inversion (findsv): "-N^ins5<invMID>ins3[EXTRA]". ToVarsBuilder renders it
                // as <INV> with genotype2 "<INV{N}>". The INV_NUM / SOME_SV_NUMBERS matches drive the
                // rendering: refallele = single base at position (no 5' anchor since '^' is present),
                // startPosition = position, endPosition = position + N - 1. MSI as proceedVrefIsDeletion.
                int dl = std::stoi(allele.substr(1));
                var.varallele = "<INV>";
                var.refallele = std::string(1, refBase);
                var.startPosition = position;
                var.endPosition = position + dl - 1;
                if (dl < cfg.SVMINLEN) {
                    // proceedVrefIsDeletion MSI (same as the plain sub-SV deletion path below).
                    std::string leftseq, tseq;
                    for (int q = std::max(position - 70, 1); q <= position - 1; ++q) if (ref.has(q)) leftseq += ref.at(q);
                    for (int q = position; q <= position + dl + 70; ++q) if (ref.has(q)) tseq += ref.at(q);
                    std::string t1 = tseq.substr(0, std::min((size_t)dl, tseq.size()));
                    std::string t2 = tseq.size() > (size_t)dl ? tseq.substr(dl) : std::string();
                    MSIResult m = findMSI(t1, t2, leftseq);
                    MSIResult m2 = findMSI(leftseq, t2);
                    double msi = m.msi; int shift3 = m.shift3; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; }
                    if (dl > 0 && msi <= (double)shift3 / dl) msi = (double)shift3 / dl;
                    var.msi = msi; var.shift3 = shift3; var.msint = msint;
                }
            } else if (allele[0] == '-') {             // deletion "-N" (possibly with complex tail &/#/^)
                int dl = std::stoi(allele.substr(1));
                bool complexDel = allele.find('&') != std::string::npos ||
                                  allele.find('#') != std::string::npos ||
                                  allele.find('^') != std::string::npos;
                // VarDict (ToVarsBuilder) anchors a *simple* deletion one base 5' of `position`:
                //   varallele = ref[position-1];  refallele = ref[position-1] + ref[position..position+dl-1]
                //   startPosition-- ; endPosition = position + dl - 1.
                char anchor = ref.has(position - 1) ? ref.at(position - 1) : refBase;
                var.startPosition = position - 1;
                var.endPosition = position + dl - 1;
                if (dl >= cfg.SVMINLEN) {
                    // Structural deletion (deletionLength >= SVMINLEN): spelled "<DEL>". The <DEL>
                    // special case (ToVarsBuilder 796-805) sets refallele to the single base at
                    // startPosition and recomputes depth/frequency; no MSI (proceedVrefIsDeletion is
                    // not called for SVs, so shift3/MSI stay 0).
                    var.varallele = "<DEL>";
                    var.refallele = ref.has(var.startPosition) ? std::string(1, ref.at(var.startPosition)) : "";
                    int tpc = var.totalPosCoverage;
                    auto cprev = vd.refCoverage.find(var.startPosition - 1);
                    if (cprev != vd.refCoverage.end()) tpc = cprev->second;
                    if (v.varsCount > tpc) tpc = v.varsCount;
                    var.totalPosCoverage = tpc;
                    var.frequency = tpc > 0 ? (double)v.varsCount / tpc : 0;
                } else {
                    // proceedVrefIsDeletion: MSI over deleted unit vs flanks (leftseq = ref[p-70..p-1],
                    // tseq = ref[p..p+dl+70]; findMSI(tseq[0..dl), tseq[dl..], leftseq) vs without-left).
                    // Java computes this for every sub-SV deletion, complex or not.
                    std::string leftseq, tseq;
                    for (int q = std::max(position - 70, 1); q <= position - 1; ++q) if (ref.has(q)) leftseq += ref.at(q);
                    for (int q = position; q <= position + dl + 70; ++q) if (ref.has(q)) tseq += ref.at(q);
                    std::string t1 = tseq.substr(0, std::min((size_t)dl, tseq.size()));
                    std::string t2 = tseq.size() > (size_t)dl ? tseq.substr(dl) : std::string();
                    MSIResult m = findMSI(t1, t2, leftseq);
                    MSIResult m2 = findMSI(leftseq, t2);
                    double msi = m.msi; int shift3 = m.shift3; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; } // shift3 unchanged
                    if (dl > 0 && msi <= (double)shift3 / dl) msi = (double)shift3 / dl;
                    var.msi = msi; var.shift3 = shift3; var.msint = msint;
                    if (!complexDel) {
                        std::string delseq;
                        for (int i = 0; i < dl; ++i) if (ref.has(position + i)) delseq += ref.at(position + i);
                        var.varallele = std::string(1, anchor);
                        var.refallele = std::string(1, anchor) + delseq;
                    } else {
                        // Complex deletion "-N&ss"/"-N^ins"/"-N#seg^tail": the 5' anchor is NOT
                        // prepended (Java skips the startPosition-- block when the description
                        // contains &/#/^). Seed refallele with the deleted bases, varallele with the
                        // tail, then apply the AMP_ATGC/HASH_CARET/CARET grammar (ToVarsBuilder 765-846).
                        std::string varAll = allele.substr(1);       // drop '-'
                        { size_t z = 0; while (z < varAll.size() && isdigit((unsigned char)varAll[z])) z++; varAll = varAll.substr(z); }
                        std::string refAll;
                        for (int i = position; i <= position + dl - 1; ++i) if (ref.has(i)) refAll += ref.at(i);
                        int sp = position, ep = position + dl - 1;
                        std::string dummy;
                        applyComplexGrammar(allele, ref, refAll, varAll, sp, ep, dummy);
                        var.refallele = refAll; var.varallele = varAll;
                        var.startPosition = sp; var.endPosition = ep;
                    }
                }
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
            // SV_info column ("<splits>-<pairs>-<clusters>", ToVarsBuilder line 428).
            if (isSVpos) var.svInfo = std::to_string(svIt->second.splits) + "-" +
                                      std::to_string(svIt->second.pairs) + "-" +
                                      std::to_string(svIt->second.clusters);
            // VarDict genotype = genotype1 + "/" + genotype2 (ToVarsBuilder). genotype1 is the
            // *description string* of the reference base (when ref freq >= -f) else the dominant
            // variant; genotype2 is THIS allele's raw description: SNV -> base, deletion -> "-N",
            // MNV -> allele with '&' stripped. Then '&'/'#' removed and '^' -> 'i'.
            {
                // genotype2 = THIS allele's raw description: deletion "-N", insertion "+len",
                // MNV/SNV keep the raw string (incl. '&', cleaned at the very end).
                auto rawDesc = [&](const std::string& a) -> std::string {
                    if (!a.empty() && a[0] == '+') return "+" + std::to_string((int)a.size() - 1);
                    if (a.size() > 1 && a[0] == '-' && a.find("<inv") != std::string::npos)
                        return "<INV" + a.substr(1, a.find('^') - 1) + ">"; // split-read INV -> "<INV{N}>"
                    return a; // deletion "-N" or MNV/SNV description (with '&' if present)
                };
                // genotype1current starts from the position's dominant description string.
                std::string g1 = positionGenotype1;
                std::string g2 = rawDesc(allele);

                // AMP_ATGC ("&([ATGC]+)"): when the variant description is followed by a matched
                // sequence, append the corresponding reference bases to genotype1current and refallele
                // and advance endPosition (ToVarsBuilder lines 765-787). The local endPosition here is
                // `position` for SNP/MNP, or position+dellen-1 for a deletion description - NOT the
                // realized var.endPosition set above for MNVs.
                int gEnd = position;
                if (!allele.empty() && allele[0] == '-') {
                    int dl = 0; size_t k = 1; while (k < allele.size() && isdigit((unsigned char)allele[k])) { dl = dl*10 + (allele[k]-'0'); ++k; }
                    gEnd = position + dl - 1;
                }
                auto ampAppend = [&](const std::string& s) -> std::string {
                    // return the matched group of the first "&([ATGC]+)" in s, else empty
                    auto amp = s.find('&');
                    if (amp == std::string::npos) return std::string();
                    std::string grp;
                    for (size_t i = amp + 1; i < s.size(); ++i) {
                        char c = s[i];
                        if (c=='A'||c=='T'||c=='G'||c=='C') grp += c; else break;
                    }
                    return grp;
                };
                auto joinRefLocal = [&](int from, int to) -> std::string {
                    std::string r;
                    for (int i = from; i <= to; ++i) if (ref.has(i)) r += ref.at(i);
                    return r;
                };
                std::string extra = ampAppend(allele);
                if (!extra.empty()) {
                    std::string tch = joinRefLocal(gEnd + 1, gEnd + (int)extra.size());
                    g1 += tch;
                    gEnd += (int)extra.size();
                    // nested AMP_ATGC on the (single-'&'-stripped) variant allele
                    std::string va = allele; auto amp = va.find('&'); if (amp != std::string::npos) va.erase(amp, 1);
                    std::string vextra = ampAppend(va);
                    if (!vextra.empty()) {
                        std::string tch2 = joinRefLocal(gEnd + 1, gEnd + (int)vextra.size());
                        g1 += tch2;
                        gEnd += (int)vextra.size();
                    }
                }

                std::string genotype = g1 + "/" + g2;
                std::string cleaned;
                for (char ch : genotype) { if (ch == '&' || ch == '#') continue; cleaned += (ch == '^') ? 'i' : ch; }
                var.genotype = cleaned;
            }

            // Reference-context flanks: 20 bp windows (ToVarsBuilder REF_20_BASES).
            for (int i = 20; i >= 1; --i) if (var.startPosition - i >= 1 && ref.has(var.startPosition - i)) var.leftseq += ref.at(var.startPosition - i);
            for (int i = 1; i <= 20; ++i) if (ref.has(var.endPosition + i)) var.rightseq += ref.at(var.endPosition + i);

            // isGoodVar sees the Java-rounded frequency (Variant.isGoodVar tests the *formatted* AF):
            // a variant at exactly -f (raw 0.009966 -> 0.0100) passes the `frequency < -f` gate, matching
            // VarDictJava, while the stored var.frequency stays raw for identical output formatting.
            { Variant probe = var; roundProbeForFilter(probe);
              if (!cfg.doPileup && !isGoodVar(cfg, probe, refHicnt, roundN(refMeanMapq, 1), vd.splice)) continue; }
            if (var.vartype == "Complex") adjComplexVar(var);   // SimplePostProcessModule.adjComplex
            result.push_back(std::move(var));
        }

        // Insertions anchored at this position.
        auto insIt = vd.insertionVariants.find(position);
        if (insIt != vd.insertionVariants.end()) {
            // ToVarsBuilder.createInsertion mutates a running hicov across the sorted insertion
            // loop: `if (hicov < hicnt) hicov = hicnt;`. Seeded with the position's calcHicov value,
            // it is a monotonic max over the insertion hi-qual counts, so an insertion whose hi-qual
            // read count exceeds the non-insertion hi-qual coverage gets hicov = hicnt (HiAF = 1.0).
            int insHicov = hicov;
            for (const auto& [allele, v] : insIt->second) {   // std::map sorted == Java Collections.sort
                // createInsertion bumps hicov for EVERY insertion (no minReads skip there); the
                // minReads/isGoodVar drop happens later, so the running max must see all insertions.
                if (insHicov < v.highQualityReadsCount) insHicov = v.highQualityReadsCount;
                if (v.varsCount < cfg.minReads) continue;
                // ToVarsBuilder.createInsertion divides by the per-insertion ttcov, NOT the position
                // Depth: when varsCount > totalCov (over-coverage) ttcov rises to varsCount so AF caps at
                // 1.0 instead of exceeding it.
                auto ttIt = insTtcov.find(allele);
                int insCov = (ttIt != insTtcov.end()) ? ttIt->second : totalCov;
                double af = insCov > 0 ? (double)v.varsCount / (double)insCov : 0.0;
                // Per-variant -f filter removed: the position-level maxfreq gate above already
                // decided whether this position survives (ToVarsBuilder emits every insertion at a
                // surviving position, subject to isGoodVar).
                Variant var;
                var.startPosition = position; var.endPosition = position;
                var.totalPosCoverage = finalTotalCov; var.varsCount = v.varsCount;
                var.varFwd = v.varsCountOnForward; var.varRev = v.varsCountOnReverse;
                var.refFwd = refFwdOut; var.refRev = refRevOut;
                var.frequency = af;
                // ToVarsBuilder.createInsertion: extraFrequency = extracnt / ttcov (same denominator as
                // frequency). realignins accumulates extracnt into the insertion Variation via adjCnt;
                // the insertion output loop previously never propagated it (AdjAF stuck at 0).
                var.extrafreq = (v.extracnt != 0 && insCov > 0) ? (double)v.extracnt / insCov : 0;
                var.pmean = v.varsCount ? v.meanPosition / v.varsCount : 0;
                var.qmean = v.varsCount ? v.meanQuality / v.varsCount : 0;
                var.mapq  = v.varsCount ? v.meanMappingQuality / v.varsCount : 0;
                var.nm    = v.varsCount ? v.numberOfMismatches / v.varsCount : 0;
                var.pstd  = v.pstd ? 1 : 0;
                var.qstd  = v.qstd ? 1 : 0;
                var.hicnt = v.highQualityReadsCount; var.hicov = insHicov;
                // A complex insertion ("+X&Y", "+...#...^...") carries a matched-sequence / indel tail
                // from CigarParser's M+I bridging; decode it into the realized ref/var alleles and
                // adjusted positions (ToVarsBuilder). Simple insertions keep the fast path + MSI.
                bool complexIns = allele.find('&') != std::string::npos ||
                                  allele.find('#') != std::string::npos ||
                                  allele.find("<dup") != std::string::npos;
                std::string refAll = std::string(1, refBase);
                std::string varAll = std::string(1, refBase) + allele.substr(1);
                int startPos = position, endPos = position;
                double refFreq = (refVar && totalCov > 0) ? (double)refVar->varsCount / totalCov : 0;
                std::string g2 = "+" + std::to_string((int)allele.size() - 1);
                std::string g1 = (refFreq >= freq) ? std::string(1, refBase) : g2;
                if (complexIns) {
                    applyComplexGrammar(allele, ref, refAll, varAll, startPos, endPos, g1);
                    var.vartype = classifyType(refAll, varAll);
                } else {
                    var.vartype = "Insertion";
                }
                // SV_info: an insertion anchored at an SV-marked position (e.g. realignlgins DUP)
                // shares the position-level "<splits>-<pairs>-<clusters>" string (ToVarsBuilder).
                if (isSVpos) var.svInfo = std::to_string(svIt->second.splits) + "-" +
                                          std::to_string(svIt->second.pairs) + "-" +
                                          std::to_string(svIt->second.clusters);
                var.refallele = refAll;
                var.varallele = varAll;
                var.startPosition = startPos;
                var.endPosition = endPos;
                {
                    std::string genotype = g1 + "/" + g2;
                    std::string cleaned;
                    for (char ch : genotype) { if (ch == '&' || ch == '#') continue; cleaned += (ch == '^') ? 'i' : ch; }
                    var.genotype = cleaned;
                }
                var.hifreq = insHicov > 0 ? (double)v.highQualityReadsCount / insHicov : 0;
                var.duprate = vd.duprate();
                if (!complexIns) {   // MSI is skipped for complex insertions (proceedVrefIsInsertion not called)
                    std::string ins = allele.substr(1);
                    std::string leftseq, tseq2;
                    for (int q = position - 50; q <= position; ++q) if (q >= 1) leftseq += ref.at(q);
                    for (int q = position + 1; q <= position + 70; ++q) tseq2 += ref.at(q);
                    MSIResult m = findMSI(ins, tseq2, leftseq);
                    MSIResult m2 = findMSI(leftseq, tseq2);
                    double msi = m.msi; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; }
                    if (!ins.empty() && msi <= m.shift3 / (double)ins.size()) msi = m.shift3 / (double)ins.size();
                    var.msi = msi; var.shift3 = m.shift3; var.msint = msint;
                }
                var.qratio = v.lowQualityReadsCount > 0
                           ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                           : (double)v.highQualityReadsCount / 0.5;
                var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads, cfg.bias))
                         + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads, cfg.bias));
                for (int i = 20; i >= 1; --i) if (var.startPosition - i >= 1 && ref.has(var.startPosition - i)) var.leftseq += ref.at(var.startPosition - i);
                for (int i = 1; i <= 20; ++i) if (ref.has(var.endPosition + i)) var.rightseq += ref.at(var.endPosition + i);
                { Variant probe = var; roundProbeForFilter(probe);
                  if (!cfg.doPileup && !isGoodVar(cfg, probe, refHicnt, roundN(refMeanMapq, 1), vd.splice)) continue; }
                if (var.vartype == "Complex") adjComplexVar(var);   // SimplePostProcessModule.adjComplex
                result.push_back(std::move(var));
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------------------------------
// Somatic (paired) mode variant builder. Mirrors callVariants' per-variant field computation exactly,
// but follows ToVarsBuilder in somatic/hasBam2 mode: the position-level `maxfreq <= -f` drop and the
// isGoodVar drop are NOT applied - every candidate is retained. Each variant is flagged `good`
// (isGoodVar with the Java-rounded frequency, which is why boundary variants at exactly -f survive)
// and tagged with its description string; positions are returned sorted, variants within a position
// ordered by ToVarsBuilder.sortVariants (rounded meanQuality * positionCoverage desc, desc string asc).
// The field logic is duplicated from callVariants deliberately so simple mode stays byte-for-byte
// unchanged.
std::vector<SomaticPosition> callVariantsSomatic(const Config& cfg, const Region& region,
                                                 const VariationData& vd, Reference& ref) {
    std::vector<SomaticPosition> out;
    const double freq = cfg.freq;

    std::vector<int> sortedPositions;
    sortedPositions.reserve(vd.nonInsertionVariants.size());
    for (const auto& kv : vd.nonInsertionVariants) sortedPositions.push_back(kv.first);
    std::sort(sortedPositions.begin(), sortedPositions.end());
    for (int position : sortedPositions) {
        const VarMap& alleleMap = vd.nonInsertionVariants.at(position);
        auto svIt = vd.svInfoAt.find(position);
        bool isSVpos = svIt != vd.svInfoAt.end();
        if (!isSVpos && (position < region.start || position > region.end)) continue;
        auto covIt = vd.refCoverage.find(position);
        if (covIt == vd.refCoverage.end() || covIt->second == 0) continue;
        int totalCov = covIt->second;
        char refBase = ref.at(position);

        const Variation* refVar = nullptr;
        auto rit = alleleMap.find(std::string(1, refBase));
        if (rit != alleleMap.end()) refVar = &rit->second;
        int refHicnt = refVar ? refVar->highQualityReadsCount : 0;
        double refMeanMapq = (refVar && refVar->varsCount) ? refVar->meanMappingQuality / refVar->varsCount : 0;

        int hicov = 0;
        for (const auto& [al, v] : alleleMap) hicov += v.highQualityReadsCount;

        // genotype1: identical to callVariants.
        std::string positionGenotype1;
        {
            double refFreqG = (refVar && totalCov > 0) ? (double)refVar->varsCount / totalCov : 0;
            if (refVar && refFreqG >= freq) {
                positionGenotype1 = std::string(1, refBase);
            } else {
                const std::string* best = nullptr;
                double bestKey = 0;
                auto consider = [&](const std::string& al, const Variation& vv) {
                    if (al.size() == 1 && al[0] == refBase) return;
                    if (vv.varsCount == 0) return;
                    double q = roundHalfEven("0.0", vv.meanQuality / (double)vv.varsCount);
                    double key = q * (double)vv.varsCount;
                    if (best == nullptr || key > bestKey || (key == bestKey && al < *best)) {
                        best = &al; bestKey = key;
                    }
                };
                for (const auto& [al, vv] : alleleMap) consider(al, vv);
                auto insG = vd.insertionVariants.find(position);
                if (insG != vd.insertionVariants.end())
                    for (const auto& [al, vv] : insG->second) consider(al, vv);
                if (best) positionGenotype1 = *best;
                else positionGenotype1 = std::string(1, refBase);
            }
            if (!positionGenotype1.empty() && positionGenotype1[0] == '+' &&
                positionGenotype1.find("<dup") == std::string::npos) {
                positionGenotype1 = "+" + std::to_string((int)positionGenotype1.size() - 1);
            }
        }

        SomaticPosition sp;
        sp.position = position;
        sp.refHicnt = refHicnt;
        sp.refMeanMapq = refMeanMapq;
        if (isSVpos) sp.sv = std::to_string(svIt->second.splits) + "-" +
                             std::to_string(svIt->second.pairs) + "-" +
                             std::to_string(svIt->second.clusters);
        // Reference-allele variant (Vars.referenceVariant): the ref base treated as a Variant, with the
        // same field mapping as a normal SNV variant. Consumed by the somatic LOH/StrongLOH paths.
        if (refVar) {
            Variant rv;
            rv.startPosition = position; rv.endPosition = position;
            rv.refallele = std::string(1, refBase);
            rv.varallele = std::string(1, refBase);
            rv.totalPosCoverage = totalCov;
            // The reference variant carries NO alt: its position/alt counts are 0 and the ref reads live
            // in refFwd/refRev. Its mean stats still come from the ref reads (refVar->varsCount denom).
            rv.varsCount = 0;
            rv.varFwd = 0; rv.varRev = 0;
            rv.refFwd = refVar->varsCountOnForward; rv.refRev = refVar->varsCountOnReverse;
            rv.frequency = 0;
            rv.pmean = refVar->varsCount ? refVar->meanPosition / refVar->varsCount : 0;
            rv.qmean = refVar->varsCount ? refVar->meanQuality / refVar->varsCount : 0;
            rv.mapq  = refVar->varsCount ? refVar->meanMappingQuality / refVar->varsCount : 0;
            rv.nm    = refVar->varsCount ? refVar->numberOfMismatches / refVar->varsCount : 0;
            rv.pstd  = refVar->pstd ? 1 : 0;
            rv.qstd  = refVar->qstd ? 1 : 0;
            rv.hicnt = refVar->highQualityReadsCount; rv.hicov = hicov;
            rv.hifreq = hicov > 0 ? (double)refVar->highQualityReadsCount / hicov : 0;
            rv.extrafreq = (refVar->extracnt != 0 && totalCov > 0) ? (double)refVar->extracnt / totalCov : 0;
            rv.qratio = refVar->lowQualityReadsCount > 0
                      ? (double)refVar->highQualityReadsCount / refVar->lowQualityReadsCount
                      : (double)refVar->highQualityReadsCount / 0.5;
            rv.bias = std::to_string(strandBias(rv.refFwd, rv.refRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias))
                    + ";" + std::to_string(strandBias(rv.varFwd, rv.varRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias));
            rv.duprate = vd.duprate();
            rv.vartype = classifyType(rv.refallele, rv.varallele);
            rv.descriptionString = std::string(1, refBase);
            {
                std::string genotype = positionGenotype1 + "/" + std::string(1, refBase);
                std::string cleaned;
                for (char ch : genotype) { if (ch == '&' || ch == '#') continue; cleaned += (ch == '^') ? 'i' : ch; }
                rv.genotype = cleaned;
            }
            for (int i = 20; i >= 1; --i) if (position - i >= 1 && ref.has(position - i)) rv.leftseq += ref.at(position - i);
            for (int i = 1; i <= 20; ++i) if (ref.has(position + i)) rv.rightseq += ref.at(position + i);
            sp.referenceVariant = std::move(rv);
            sp.hasRef = true;
        }

        // Records a built variant into the position group with good flag + description string.
        auto record = [&](Variant&& var, const std::string& desc) {
            var.descriptionString = desc;
            // isGoodVar sees the Java-rounded frequency so a variant at exactly -f (rounds to 0.0100)
            // is retained; the stored frequency stays raw (formatting rounds identically for output).
            Variant probe = var;
            roundProbeForFilter(probe);
            var.good = isGoodVar(cfg, probe, refHicnt, roundN(refMeanMapq, 1), vd.splice);
            sp.variants.push_back(std::move(var));
        };

        for (const auto& [allele, v] : alleleMap) {
            if (allele.size() == 1 && allele[0] == refBase) continue;
            // Somatic keeps every allele (even below -r) so the cross-sample descriptionString lookups
            // (LOH / paired match) can find it; -r only gates the good flag via isGoodVar.
            double af = totalCov > 0 ? (double)v.varsCount / (double)totalCov : 0.0;

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
            var.extrafreq = (v.extracnt != 0 && totalCov > 0) ? (double)v.extracnt / totalCov : 0;
            var.qratio = v.lowQualityReadsCount > 0
                       ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                       : (double)v.highQualityReadsCount / 0.5;
            var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias))
                     + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads == 0 ? 2 : cfg.minBiasReads, cfg.bias));
            var.duprate = vd.duprate();

            if (allele.find('&') != std::string::npos && allele[0] != '-' && allele[0] != '+') {   // MNV / complex (e.g. "A&CG" -> ACG)
                std::string va;
                for (char ch : allele) if (ch != '&') va += ch;
                var.varallele = va;
                std::string ra;
                for (int k = 0; k < (int)va.size(); ++k) ra += ref.at(position + k);
                var.refallele = ra;
                var.endPosition = position + (int)va.size() - 1;
                std::string t1, t2;
                for (int q = position - 30; q <= position + 1; ++q) if (q >= 1) t1 += ref.at(q);
                for (int q = position + 2; q <= position + 70; ++q) t2 += ref.at(q);
                MSIResult m = findMSI(t1, t2);
                var.msi = m.msi; var.shift3 = m.shift3; var.msint = m.msintLen;
            } else if (allele[0] == '+') {
                var.refallele = std::string(1, refBase);
                var.varallele = std::string(1, refBase) + allele.substr(1);
            } else if (allele[0] == '-' && allele.find("<inv") != std::string::npos) {
                // Split-read inversion (findsv): rendered as <INV> with genotype2 "<INV{N}>".
                int dl = std::stoi(allele.substr(1));
                var.varallele = "<INV>";
                var.refallele = std::string(1, refBase);
                var.startPosition = position;
                var.endPosition = position + dl - 1;
                if (dl < cfg.SVMINLEN) {
                    std::string leftseq, tseq;
                    for (int q = std::max(position - 70, 1); q <= position - 1; ++q) if (ref.has(q)) leftseq += ref.at(q);
                    for (int q = position; q <= position + dl + 70; ++q) if (ref.has(q)) tseq += ref.at(q);
                    std::string t1 = tseq.substr(0, std::min((size_t)dl, tseq.size()));
                    std::string t2 = tseq.size() > (size_t)dl ? tseq.substr(dl) : std::string();
                    MSIResult m = findMSI(t1, t2, leftseq);
                    MSIResult m2 = findMSI(leftseq, t2);
                    double msi = m.msi; int shift3 = m.shift3; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; }
                    if (dl > 0 && msi <= (double)shift3 / dl) msi = (double)shift3 / dl;
                    var.msi = msi; var.shift3 = shift3; var.msint = msint;
                }
            } else if (allele[0] == '-') {             // deletion "-N" (possibly with complex tail &/#/^)
                int dl = std::stoi(allele.substr(1));
                bool complexDel = allele.find('&') != std::string::npos ||
                                  allele.find('#') != std::string::npos ||
                                  allele.find('^') != std::string::npos;
                char anchor = ref.has(position - 1) ? ref.at(position - 1) : refBase;
                var.startPosition = position - 1;
                var.endPosition = position + dl - 1;
                if (dl >= cfg.SVMINLEN) {
                    var.varallele = "<DEL>";
                    var.refallele = ref.has(var.startPosition) ? std::string(1, ref.at(var.startPosition)) : "";
                    int tpc = var.totalPosCoverage;
                    auto cprev = vd.refCoverage.find(var.startPosition - 1);
                    if (cprev != vd.refCoverage.end()) tpc = cprev->second;
                    if (v.varsCount > tpc) tpc = v.varsCount;
                    var.totalPosCoverage = tpc;
                    var.frequency = tpc > 0 ? (double)v.varsCount / tpc : 0;
                } else {
                    std::string leftseq, tseq;
                    for (int q = std::max(position - 70, 1); q <= position - 1; ++q) if (ref.has(q)) leftseq += ref.at(q);
                    for (int q = position; q <= position + dl + 70; ++q) if (ref.has(q)) tseq += ref.at(q);
                    std::string t1 = tseq.substr(0, std::min((size_t)dl, tseq.size()));
                    std::string t2 = tseq.size() > (size_t)dl ? tseq.substr(dl) : std::string();
                    MSIResult m = findMSI(t1, t2, leftseq);
                    MSIResult m2 = findMSI(leftseq, t2);
                    double msi = m.msi; int shift3 = m.shift3; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; }
                    if (dl > 0 && msi <= (double)shift3 / dl) msi = (double)shift3 / dl;
                    var.msi = msi; var.shift3 = shift3; var.msint = msint;
                    if (!complexDel) {
                        std::string delseq;
                        for (int i = 0; i < dl; ++i) if (ref.has(position + i)) delseq += ref.at(position + i);
                        var.varallele = std::string(1, anchor);
                        var.refallele = std::string(1, anchor) + delseq;
                    } else {
                        // Complex deletion "-N&ss"/"-N^ins"/"-N#seg^tail": the 5' anchor is NOT
                        // prepended; seed refallele with deleted bases, varallele with the tail, then
                        // apply the AMP_ATGC/HASH_CARET/CARET grammar (ToVarsBuilder 765-846).
                        std::string varAll = allele.substr(1);       // drop '-'
                        { size_t z = 0; while (z < varAll.size() && isdigit((unsigned char)varAll[z])) z++; varAll = varAll.substr(z); }
                        std::string refAll;
                        for (int i = position; i <= position + dl - 1; ++i) if (ref.has(i)) refAll += ref.at(i);
                        int sp2 = position, ep = position + dl - 1;
                        std::string dummy;
                        applyComplexGrammar(allele, ref, refAll, varAll, sp2, ep, dummy);
                        var.refallele = refAll; var.varallele = varAll;
                        var.startPosition = sp2; var.endPosition = ep;
                    }
                }
            } else {
                var.refallele = std::string(1, refBase);
                var.varallele = allele;
                std::string tseq1, tseq2;
                for (int q = position - 30; q <= position + 1; ++q) if (q >= 1) tseq1 += ref.at(q);
                for (int q = position + 2; q <= position + 70; ++q) tseq2 += ref.at(q);
                MSIResult m = findMSI(tseq1, tseq2);
                var.msi = m.msi; var.shift3 = m.shift3; var.msint = m.msintLen;
            }
            var.vartype = classifyType(var.refallele, var.varallele);
            {
                auto rawDesc = [&](const std::string& a) -> std::string {
                    if (!a.empty() && a[0] == '+') return "+" + std::to_string((int)a.size() - 1);
                    return a;
                };
                std::string g1 = positionGenotype1;
                std::string g2 = rawDesc(allele);
                int gEnd = position;
                if (!allele.empty() && allele[0] == '-') {
                    int dl = 0; size_t k = 1; while (k < allele.size() && isdigit((unsigned char)allele[k])) { dl = dl*10 + (allele[k]-'0'); ++k; }
                    gEnd = position + dl - 1;
                }
                auto ampAppend = [&](const std::string& s) -> std::string {
                    auto amp = s.find('&');
                    if (amp == std::string::npos) return std::string();
                    std::string grp;
                    for (size_t i = amp + 1; i < s.size(); ++i) {
                        char cc = s[i];
                        if (cc=='A'||cc=='T'||cc=='G'||cc=='C') grp += cc; else break;
                    }
                    return grp;
                };
                auto joinRefLocal = [&](int from, int to) -> std::string {
                    std::string r;
                    for (int i = from; i <= to; ++i) if (ref.has(i)) r += ref.at(i);
                    return r;
                };
                std::string extra = ampAppend(allele);
                if (!extra.empty()) {
                    std::string tch = joinRefLocal(gEnd + 1, gEnd + (int)extra.size());
                    g1 += tch;
                    gEnd += (int)extra.size();
                    std::string va = allele; auto amp = va.find('&'); if (amp != std::string::npos) va.erase(amp, 1);
                    std::string vextra = ampAppend(va);
                    if (!vextra.empty()) {
                        std::string tch2 = joinRefLocal(gEnd + 1, gEnd + (int)vextra.size());
                        g1 += tch2;
                        gEnd += (int)vextra.size();
                    }
                }
                std::string genotype = g1 + "/" + g2;
                std::string cleaned;
                for (char ch : genotype) { if (ch == '&' || ch == '#') continue; cleaned += (ch == '^') ? 'i' : ch; }
                var.genotype = cleaned;
            }
            for (int i = 20; i >= 1; --i) if (var.startPosition - i >= 1 && ref.has(var.startPosition - i)) var.leftseq += ref.at(var.startPosition - i);
            for (int i = 1; i <= 20; ++i) if (ref.has(var.endPosition + i)) var.rightseq += ref.at(var.endPosition + i);

            record(std::move(var), allele);
        }

        auto insIt = vd.insertionVariants.find(position);
        if (insIt != vd.insertionVariants.end()) {
            // Same running-hicov bump as callVariants / ToVarsBuilder.createInsertion.
            int insHicov = hicov;
            for (const auto& [allele, v] : insIt->second) {   // std::map sorted == Java Collections.sort
                if (insHicov < v.highQualityReadsCount) insHicov = v.highQualityReadsCount;
                double af = totalCov > 0 ? (double)v.varsCount / (double)totalCov : 0.0;
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
                var.pstd  = v.pstd ? 1 : 0;
                var.qstd  = v.qstd ? 1 : 0;
                var.hicnt = v.highQualityReadsCount; var.hicov = insHicov;
                // A complex insertion ("+X&Y", "+...#...^...") carries a matched-sequence / indel tail
                // from CigarParser's M+I bridging; decode it into the realized ref/var alleles and
                // adjusted positions (ToVarsBuilder). Simple insertions keep the fast path + MSI.
                bool complexIns = allele.find('&') != std::string::npos ||
                                  allele.find('#') != std::string::npos ||
                                  allele.find("<dup") != std::string::npos;
                std::string refAll = std::string(1, refBase);
                std::string varAll = std::string(1, refBase) + allele.substr(1);
                int startPos = position, endPos = position;
                double refFreq = (refVar && totalCov > 0) ? (double)refVar->varsCount / totalCov : 0;
                std::string g2 = "+" + std::to_string((int)allele.size() - 1);
                std::string g1 = (refFreq >= freq) ? std::string(1, refBase) : g2;
                if (complexIns) {
                    applyComplexGrammar(allele, ref, refAll, varAll, startPos, endPos, g1);
                    var.vartype = classifyType(refAll, varAll);
                } else {
                    var.vartype = "Insertion";
                }
                // SV_info: an insertion anchored at an SV-marked position (e.g. realignlgins DUP)
                // shares the position-level "<splits>-<pairs>-<clusters>" string (ToVarsBuilder).
                if (isSVpos) var.svInfo = std::to_string(svIt->second.splits) + "-" +
                                          std::to_string(svIt->second.pairs) + "-" +
                                          std::to_string(svIt->second.clusters);
                var.refallele = refAll;
                var.varallele = varAll;
                var.startPosition = startPos;
                var.endPosition = endPos;
                {
                    std::string genotype = g1 + "/" + g2;
                    std::string cleaned;
                    for (char ch : genotype) { if (ch == '&' || ch == '#') continue; cleaned += (ch == '^') ? 'i' : ch; }
                    var.genotype = cleaned;
                }
                var.hifreq = insHicov > 0 ? (double)v.highQualityReadsCount / insHicov : 0;
                var.duprate = vd.duprate();
                if (!complexIns) {   // MSI is skipped for complex insertions (proceedVrefIsInsertion not called)
                    std::string ins = allele.substr(1);
                    std::string leftseq, tseq2;
                    for (int q = position - 50; q <= position; ++q) if (q >= 1) leftseq += ref.at(q);
                    for (int q = position + 1; q <= position + 70; ++q) tseq2 += ref.at(q);
                    MSIResult m = findMSI(ins, tseq2, leftseq);
                    MSIResult m2 = findMSI(leftseq, tseq2);
                    double msi = m.msi; int msint = m.msintLen;
                    if (msi < m2.msi) { msi = m2.msi; msint = m2.msintLen; }
                    if (!ins.empty() && msi <= m.shift3 / (double)ins.size()) msi = m.shift3 / (double)ins.size();
                    var.msi = msi; var.shift3 = m.shift3; var.msint = msint;
                }
                var.qratio = v.lowQualityReadsCount > 0
                           ? (double)v.highQualityReadsCount / v.lowQualityReadsCount
                           : (double)v.highQualityReadsCount / 0.5;
                var.bias = std::to_string(strandBias(var.refFwd, var.refRev, cfg.minBiasReads, cfg.bias))
                         + ";" + std::to_string(strandBias(var.varFwd, var.varRev, cfg.minBiasReads, cfg.bias));
                for (int i = 20; i >= 1; --i) if (var.startPosition - i >= 1 && ref.has(var.startPosition - i)) var.leftseq += ref.at(var.startPosition - i);
                for (int i = 1; i <= 20; ++i) if (ref.has(var.endPosition + i)) var.rightseq += ref.at(var.endPosition + i);
                // Insertion description strings are stored under "+<seq>" like createInsertion's key.
                record(std::move(var), allele);
            }
        }

        // sortVariants: rounded meanQuality * positionCoverage descending, tie-break description asc.
        std::sort(sp.variants.begin(), sp.variants.end(), [](const Variant& a, const Variant& b) {
            double ka = roundHalfEven("0.0", a.qmean) * a.varsCount;
            double kb = roundHalfEven("0.0", b.qmean) * b.varsCount;
            if (ka != kb) return ka > kb;
            return a.descriptionString < b.descriptionString;
        });

        out.push_back(std::move(sp));
    }
    // Positions come from a std::map (nonInsertionVariants) and are already ascending; keep it explicit.
    std::sort(out.begin(), out.end(), [](const SomaticPosition& a, const SomaticPosition& b) {
        return a.position < b.position;
    });
    return out;
}

} // namespace vardict
