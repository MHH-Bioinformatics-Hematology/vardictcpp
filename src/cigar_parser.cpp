#include "cigar_parser.hpp"
#include "util.hpp"
#include "cigar_modifier.hpp"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <string>

namespace vardict {

static inline char baseChar(const bam1_t* b, int qpos) {
    static const char code[] = "=ACMGRSVTWYHKDBN";
    return code[bam_seqi(bam_get_seq(b), qpos)];
}

// VariationUtils.getVariationFromSeq: get-or-create the Variation for a soft-clip consensus base.
static inline Variation& getVariationFromSeq(Sclip& sc, int idx, char ch) {
    auto& slot = sc.seq[idx][ch];
    if (!slot) slot = std::make_shared<Variation>();
    return *slot;
}

// VariationRealigner.adjInsPos: left-normalize an insertion anchor. Shifts `bi` left while the
// reference base equals the correspondingly-rotated insertion base, rotating `ins` accordingly, so
// alignment-ambiguous insertions in repeats collapse to one canonical (position, sequence).
static void adjInsPos(int& bi, std::string& ins, Reference& ref) {
    int n = 1;
    int len = (int)ins.size();
    while (ref.has(bi) && ref.at(bi) == ins[ins.size() - n]) {
        n++;
        if (n > len) n = 1;
        bi--;
    }
    if (n > 1) ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
}

// VariationUtils.addCnt: accumulate one observation into a Variation (no pstd/qstd; those are set
// only in the matching-part increment).
static inline void addCnt(Variation& v, bool dir, int readPosition, double baseQuality,
                          int mappingQuality, double nm, double goodq) {
    v.varsCount++;
    v.incDir(dir);
    v.meanPosition += readPosition;
    v.meanQuality += baseQuality;
    v.meanMappingQuality += mappingQuality;
    v.numberOfMismatches += nm;
    if (baseQuality >= goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
}

// CigarParser.addSV: fold one discordant read pair into an SV (deletion) cluster and record its Mate.
static void addSVMate(Sclip& sd, int start_s, int end_e, int mateStart_ms, int mateEnd_me,
                      int dir, int rlen, int mlen, int softp, double pmean_rp, double qmean,
                      double Qmean, double nm, double goodq) {
    sd.varsCount++;
    sd.incDir(dir == 1 ? false : true);
    if (qmean >= goodq) sd.highQualityReadsCount++; else sd.lowQualityReadsCount++;
    if (sd.start == 0 || sd.start >= start_s) sd.start = start_s;
    if (sd.end == 0 || sd.end <= end_e) sd.end = end_e;
    Mate m; m.mateStart_ms = mateStart_ms; m.mateEnd_me = mateEnd_me; m.mateLength_mlen = mlen;
    m.start_s = start_s; m.end_e = end_e; m.pmean_rp = pmean_rp; m.qmean_q = qmean; m.Qmean_Q = Qmean; m.nm = nm;
    sd.mates.push_back(m);
    if (sd.mstart == 0 || sd.mstart >= mateStart_ms) sd.mstart = mateStart_ms;
    if (sd.mend == 0 || sd.mend <= mateEnd_me) sd.mend = mateStart_ms + rlen; // NB: mateStart+rlen, not mateEnd_me (VarDict addSV:2322)
    if (softp != 0) {
        if (dir == 1) { if (std::abs(softp - sd.end) < 10) sd.soft[softp]++; }
        else          { if (std::abs(softp - sd.start) < 10) sd.soft[softp]++; }
    }
}

// CigarParser.prepareSVStructuresForAnalysis, DELETION (svfdel/svrdel) path only. Records a
// same-chromosome discordant read pair whose orientation + oversized insert signal a deletion into the
// forward/reverse DEL cluster list, opening a new cluster when the read is > MINSVCDIST*maxReadLength
// past the last one. DUP/INV/inter-chromosome classification is not ported (no golden), so only the DEL
// disc-count bumps are applied; that can only lower a DEL cluster's disc, never change its membership.
static void prepareSVDel(const bam1_t* b, const Cig& cigv, int start,
                         const std::vector<int>& bqual, int lqseq, bool reverse,
                         double nm, VariationData& out, const Config& cfg) {
    const bam1_core_t& c = b->core;
    if (c.mtid != c.tid) return; // getMateReferenceName == "=" (same chr); inter-chr not ported
    int totalLen = 0, alnMND = 0;
    for (auto& e : cigv) {
        char op = e.second;
        if (op=='M'||op=='I'||op=='S'||op=='='||op=='X') totalLen += e.first;
        if (op=='M'||op=='N'||op=='D') alnMND += e.first;
    }
    int end = start + alnMND;
    int mateStart = c.mpos + 1;
    int mend = mateStart + totalLen;
    int soft5 = 0, soft3 = 0;
    if (!cigv.empty() && cigv.front().second == 'S') {
        int tt = cigv.front().first;
        if (tt != 0 && tt - 1 < lqseq && bqual[tt - 1] > cfg.goodq) soft5 = start;
    }
    if (!cigv.empty() && cigv.back().second == 'S') {
        int tt = cigv.back().first;
        if (tt != 0 && lqseq - tt >= 0 && lqseq - tt < lqseq && bqual[lqseq - tt] > cfg.goodq) soft3 = end;
    }
    int readDirNum = reverse ? -1 : 1;
    bool mateForward = (c.flag & 0x20) == 0;
    int mateDirNum = mateForward ? 1 : -1;
    long mlen = c.isize;
    if (uint8_t* mc = bam_aux_get(const_cast<bam1_t*>(b), "MC")) {
        const char* s = bam_aux2Z(mc);
        if (s) { int cntS = 0; for (const char* p = s; *p; ++p) if (*p == 'S') cntS++; if (cntS >= 2) return; }
    }
    if (uint8_t* mq = bam_aux_get(const_cast<bam1_t*>(b), "MQ")) {
        if ((int)bam_aux2i(mq) < 15) return;
    }
    if (readDirNum * mateDirNum == -1 && (mlen * readDirNum) > 0 && lqseq > Config::MINMAPBASE) {
        mlen = mateStart > start ? (long)mend - start : (long)end - mateStart;
        if (std::labs(mlen) > (long)cfg.INSSIZE + (long)cfg.INSSTDAMT * cfg.INSSTD) {
            double qAtBase = bqual[Config::MINMAPBASE];
            double Qmean = c.qual;
            double pmean = out.maxReadLength / 2.0;
            if (readDirNum == 1) {
                if (out.svfdel.empty() || start - out.svdelfend > Config::MINSVCDIST * out.maxReadLength)
                    { Sclip sc; sc.varsCount = 0; out.svfdel.push_back(sc); }
                addSVMate(out.svfdel.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft3, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svdelfend = end;
            } else {
                if (out.svrdel.empty() || start - out.svdelrend > Config::MINSVCDIST * out.maxReadLength)
                    { Sclip sc; sc.varsCount = 0; out.svrdel.push_back(sc); }
                addSVMate(out.svrdel.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft5, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svdelrend = end;
            }
            if (!out.svfdel.empty() && std::abs(start - out.svdelfend) <= Config::MINSVCDIST * out.maxReadLength)
                out.svfdel.back().disc++;
            if (!out.svrdel.empty() && std::abs(start - out.svdelrend) <= Config::MINSVCDIST * out.maxReadLength)
                out.svrdel.back().disc++;
        }
    }
}

bool CigarParser::process(const Region& region, VariationData& out) {
    samFile* fp = sam_open(cfg_.bam.c_str(), "r");
    if (!fp) throw std::runtime_error("cannot open BAM " + cfg_.bam);
    hts_idx_t* idx = sam_index_load(fp, cfg_.bam.c_str());
    if (!idx) throw std::runtime_error("cannot load BAM index for " + cfg_.bam + " (run `samtools index`)");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) throw std::runtime_error("cannot read BAM header");

    int tid = bam_name2id(hdr, region.chr.c_str());
    if (tid < 0) { // try toggling chr prefix
        std::string alt = (region.chr.rfind("chr", 0) == 0) ? region.chr.substr(3) : "chr" + region.chr;
        tid = bam_name2id(hdr, alt.c_str());
    }
    if (tid < 0) { bam_hdr_destroy(hdr); hts_idx_destroy(idx); sam_close(fp); return true; }

    out.chrLen = hdr->target_len[tid];
    // Fetch reads overlapping the region (htslib is 0-based, end-exclusive).
    hts_itr_t* it = sam_itr_queryi(idx, tid, region.start - 1, region.end);
    bam1_t* b = bam_init1();

    const int rlo = region.start, rhi = region.end;
    // Simple duplicate accounting for the duprate column (start+cigar signature per alignment start).
    long lastStart = -1;
    std::map<std::string, int> dupKeys;

    while (sam_itr_next(fp, it, b) >= 0) {
        const bam1_core_t& c = b->core;
        if (c.flag & BAM_FUNMAP) continue;
        if (cfg_.samFilterFlag != 0 && (c.flag & cfg_.samFilterFlag)) continue;
        if ((double)c.qual < cfg_.mapqMin) continue;
        if (c.n_cigar == 0 || c.l_qseq == 0) continue;
        // Ignore supplementary alignments (they would skew coverage), as VarDict does.
        if (c.flag & BAM_FSUPPLEMENTARY) continue;
        // Ignore reads soft-clipped at both ends where the leading clip is 10-99 bp and the trailing
        // clip is >= 10 bp (VarDict pattern ^\d\dS.*\d\dS$ on the CIGAR: two-digit leading S, >=2-digit
        // trailing S). These are chimeric/mis-mapped islands VarDict does not count.
        {
            const uint32_t* cg0 = bam_get_cigar(b);
            int lead = (bam_cigar_op(cg0[0]) == BAM_CSOFT_CLIP) ? (int)bam_cigar_oplen(cg0[0]) : 0;
            int tail = (bam_cigar_op(cg0[c.n_cigar - 1]) == BAM_CSOFT_CLIP) ? (int)bam_cigar_oplen(cg0[c.n_cigar - 1]) : 0;
            if (lead >= 10 && lead <= 99 && tail >= 10) continue;
        }

        bool reverse = (c.flag & BAM_FREVERSE) != 0;
        int mapq = c.qual;

        // NM tag (edit distance) if present. VarDict (CigarParser) subtracts the total I+D length
        // of the ORIGINAL cigar so per-variation NM counts mismatches only, not the indel gaps
        // ("Edit distance - indels is the # of mismatches"). Computed before modifyCigar.
        double nm = 0;
        bool haveNM = false;
        if (uint8_t* aux = bam_aux_get(b, "NM")) {
            long insDelLen = 0;
            const uint32_t* nmcig = bam_get_cigar(b);
            for (uint32_t k = 0; k < c.n_cigar; ++k) {
                int op = bam_cigar_op(nmcig[k]);
                if (op == BAM_CINS || op == BAM_CDEL) insDelLen += bam_cigar_oplen(nmcig[k]);
            }
            nm = (double)bam_aux2i(aux) - (double)insDelLen;
            haveNM = true;
        }
        // VarDict (CigarParser): "reads with mismatches more than INT will be filtered and ignored"
        // (gaps not counted). Skip the whole read when NM - indels exceeds -m (default 8).
        if (haveNM && nm > cfg_.mismatch) continue;

        // Amplicon read-assignment (CigarParser.parseCigarWithAmpCase): keep a read only if it belongs
        // to THIS amplicon (region == the amplicon interval). Applied here, before maxReadLength/counting,
        // matching Java's order (NM filter -> amp check -> modifyCigar -> maxReadLength). Uses the ORIGINAL
        // cigar. Returns skip=true to drop the read.
        if (cfg_.amplicon) {
            const int distanceToAmplicon = cfg_.ampEdge;
            const double overlapFraction = cfg_.ampFraction;
            const uint32_t* acg = bam_get_cigar(b);
            // getAlignedLength: sum of M and D lengths (aligned length, excl. soft-clip and insertion).
            int alignedLen = 0;
            for (uint32_t k = 0; k < c.n_cigar; ++k) {
                int op = bam_cigar_op(acg[k]);
                if (op == BAM_CMATCH || op == BAM_CDEL) alignedLen += (int)bam_cigar_oplen(acg[k]);
            }
            int segstart = c.pos + 1;                 // record.getAlignmentStart(), 1-based
            int segend = segstart + alignedLen - 1;
            int firstOp = bam_cigar_op(acg[0]);
            int lastOp  = bam_cigar_op(acg[c.n_cigar - 1]);
            bool skip = false;
            if (firstOp == BAM_CSOFT_CLIP) {
                int ts1 = segstart > region.start ? segstart : region.start;
                int te1 = segend < region.end ? segend : region.end;
                if (!(std::abs(ts1 - te1) / (double)(segend - segstart) > overlapFraction)) skip = true;
            } else if (lastOp == BAM_CSOFT_CLIP) {
                int ts1 = segstart > region.start ? segstart : region.start;
                int te1 = segend < region.end ? segend : region.end;
                if (!(std::abs(te1 - ts1) / (double)(segend - segstart) > overlapFraction)) skip = true;
            } else { // no soft-clipping: use mate/TLEN to bound the fragment
                bool isMateReferenceNameEqual = (c.mtid == c.tid);
                if (isMateReferenceNameEqual && c.isize != 0) {
                    if (c.isize > 0) segend = segstart + c.isize - 1;
                    else { segstart = c.mpos + 1; segend = c.mpos + 1 - c.isize - 1; }
                }
                int ts1 = segstart > region.start ? segstart : region.start;
                int te1 = segend < region.end ? segend : region.end;
                if ((std::abs(segstart - region.start) > distanceToAmplicon ||
                     std::abs(segend - region.end) > distanceToAmplicon)
                    || std::fabs((ts1 - te1) / (double)(segend - segstart)) <= overlapFraction) {
                    skip = true;
                }
            }
            if (skip) continue;
        }

        out.totalReads++;
        if (c.l_qseq > out.maxReadLength) out.maxReadLength = c.l_qseq;

        // Duplicate detection (optional): same start + cigar seen again.
        if (cfg_.removeDuplicates) {
            if (c.pos != lastStart) { dupKeys.clear(); lastStart = c.pos; }
            std::string key = std::to_string(c.pos) + ":";
            const uint32_t* cg = bam_get_cigar(b);
            for (uint32_t k = 0; k < c.n_cigar; ++k) key += std::to_string(cg[k]) + ",";
            if (++dupKeys[key] > 1) { out.dupReads++; continue; }
        }

        // Build read sequence + qualities + CIGAR, then reshape the CIGAR (CigarModifier) before counting.
        static const char CODE[] = "=ACMGRSVTWYHKDBN";
        const uint8_t* rawseq = bam_get_seq(b);
        const uint8_t* rawqual = bam_get_qual(b);
        std::string bseq(c.l_qseq, 'N');
        std::vector<int> bqual(c.l_qseq);
        for (int j = 0; j < c.l_qseq; ++j) { bseq[j] = CODE[bam_seqi(rawseq, j)]; bqual[j] = rawqual[j]; }
        const uint32_t* rawcig = bam_get_cigar(b);
        static const char OPS[] = "MIDNSHP=X";
        Cig cigv;
        for (uint32_t k = 0; k < c.n_cigar; ++k) cigv.push_back({(int)bam_cigar_oplen(rawcig[k]), OPS[bam_cigar_op(rawcig[k])]});
        int rpos = c.pos + 1;   // 1-based reference position of current op
        if (cfg_.performLocalRealignment) modifyCigar(rpos, cigv, bseq, bqual, ref_, out.maxReadLength, cfg_);

        // Structural-variant discordant-pair collection (CigarParser dispatch at 323-329): skip
        // paired reads whose mate is unmapped (potential insertion, not ported); otherwise, for
        // MAPQ>10 reads, record possible SV deletion clusters from the (post-modifyCigar) alignment.
        if (!cfg_.disableSV) {
            bool paired = (c.flag & BAM_FPAIRED) != 0;
            bool mateUnmapped = (c.flag & BAM_FMUNMAP) != 0;
            if (paired && mateUnmapped) { /* potential insertion: not ported */ }
            else if (c.qual > 10) prepareSVDel(b, cigv, rpos, bqual, c.l_qseq, reverse, nm, out, cfg_);
        }

        int qpos = 0;           // 0-based query offset (includes soft-clip)

        // VarDict read-position convention: position within the aligned read (M+I only, excluding
        // soft-clip), folded to the distance from the nearest read end. Precompute the aligned length.
        int rlen = 0;           // readLengthIncludeMatchingAndInsertions
        for (auto& e : cigv) if (e.second == 'M' || e.second == '=' || e.second == 'X' || e.second == 'I') rlen += e.first;
        int rpe = 0;            // readPositionExcludingSoftClipped
        auto foldPos = [&](int posExclSc) {
            return posExclSc < rlen - posExclSc ? posExclSc + 1 : rlen - posExclSc;
        };

        for (uint32_t k = 0; k < cigv.size(); ++k) {
            char op = cigv[k].second;
            int len = cigv[k].first;
            switch (op) {
            case 'M':
            case '=':
            case 'X': {
                // Faithful port of CigarParser's matching-part loop with MNV growth: adjacent
                // mismatches (bridging up to vext matching bases) are grown into a single variant
                // whose description string joins the leading base(s) and the grown tail with '&'
                // (e.g. "A&CG"). Single-base matches (the common case) reduce to one variation "X".
                // Adjacent-indel bridging within an M-segment (the '^'/'#'/'-N&' grammar) is not
                // grown here; those complex cases are left to the I/D handlers.
                int i = 0;
                // nmoff accumulates over the WHOLE M-segment (declared once, like CigarParser.java):
                // each consecutive-mismatch base absorbed into an MNV increments it, and every
                // variation in this segment contributes (nm - nmoff) so the read's per-base mismatches
                // that were merged into an MNV are not double-counted as "other" mismatches.
                int nmoff = 0;
                while (i < len) {
                    int gref = rpos + i;         // moving reference position (VarDict 'start')
                    int gq   = qpos + i;         // moving query offset (incl. soft-clip)
                    int grpe = rpe + i;          // moving readPositionExcludingSoftClipped
                    char ch1 = bseq[gq];
                    if (ch1 == 'N') { i++; continue; }
                    double q = bqual[gq];
                    int qbases = 1;
                    std::string s(1, ch1);
                    std::string ss;
                    // Grow MNV while the current base mismatches the reference and quality is good.
                    while ((gref + 1) >= rlo && (gref + 1) <= rhi && (i + 1) < len &&
                           q >= cfg_.goodq &&
                           ref_.at(gref) != bseq[gq] && ref_.at(gref) != 'N') {
                        if (bqual[gq + 1] < cfg_.goodq + 5) break;
                        char nuc = bseq[gq + 1];
                        if (nuc == 'N') break;
                        if (ref_.at(gref + 1) == 'N') break;
                        if (ref_.at(gref + 1) != nuc) {           // next base also mismatches
                            ss += nuc; q += bqual[gq + 1]; qbases++;
                            gq++; gref++; grpe++; i++;
                            nmoff++;                              // CigarParser.java: nmoff++ per absorbed consecutive mismatch
                        } else {                                  // bridge matching bases to next mismatch within vext
                            int ssn = 0;
                            for (int ssi = 1; ssi <= cfg_.vext; ssi++) {
                                if (i + 1 + ssi >= len) break;
                                if (bseq[gq + 1 + ssi] != ref_.at(gref + 1 + ssi)) { ssn = ssi + 1; break; }
                            }
                            if (ssn == 0) break;
                            if (bqual[gq + ssn] < cfg_.goodq + 5) break;
                            for (int ssi = 1; ssi <= ssn; ssi++) { ss += bseq[gq + ssi]; q += bqual[gq + ssi]; qbases++; }
                            gq += ssn; gref += ssn; grpe += ssn; i += ssn;
                        }
                    }
                    if (!ss.empty()) s += "&" + ss;
                    int pos = gref - qbases + 1;                  // leftmost covered position
                    double qavg = q / qbases;
                    int tp = grpe < rlen - grpe ? grpe + 1 : rlen - grpe;
                    if (pos >= rlo && pos <= rhi && s.find('N') == std::string::npos) {
                        Variation& v = out.nonInsertionVariants[pos][s];
                        if (!v.pstd && v.pp != 0 && tp != v.pp) v.pstd = true;
                        if (!v.qstd && v.pq != 0 && qavg != v.pq) v.qstd = true;
                        v.varsCount++;
                        v.incDir(reverse);
                        v.meanPosition += tp;
                        v.meanQuality += qavg;
                        v.meanMappingQuality += mapq;
                        v.numberOfMismatches += nm - nmoff;      // subtract mismatches merged into MNVs earlier in this segment
                        v.pp = tp; v.pq = qavg;
                        if (qavg >= cfg_.goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
                        // reference coverage for every base covered by this variation
                        for (int qi = 1; qi <= qbases; ++qi) {
                            int cp = gref - qi + 1;
                            if (cp >= rlo && cp <= rhi) out.refCoverage[cp]++;
                        }
                        // MNP bookkeeping (one base + '&' + more bases)
                        if (ss.size() >= 1 && s.find('&') != std::string::npos)
                            out.mnp[pos][s]++;
                    }
                    i++;
                }
                rpos += len; qpos += len; rpe += len;
                break;
            }
            case 'I': {
                int p = rpos - 1; // insertion anchored to preceding reference base (VarDict convention)
                std::string ins;
                double qsum = 0;
                for (int i = 0; i < len; ++i) { ins += bseq[qpos + i]; qsum += bqual[qpos + i]; }
                adjInsPos(p, ins, ref_);   // left-normalize insertion anchor in repeats
                std::string sig = "+" + ins;
                if (p >= rlo && p <= rhi && ins.find('N') == std::string::npos) {
                    int tp = foldPos(rpe);
                    // mean quality of the inserted segment (CigarParser.processInsertion: tmpq)
                    double tmpq = qsum / (len ? len : 1);
                    out.positionToInsertionCount[p][sig]++;
                    Variation& v = out.insertionVariants[p][sig];
                    // pstd/qstd flags (set before pp/pq are refreshed), per Java processInsertion
                    if (!v.pstd && v.pp != 0 && tp != v.pp) v.pstd = true;
                    if (!v.qstd && v.pq != 0 && tmpq != v.pq) v.qstd = true;
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += tmpq;
                    v.meanMappingQuality += mapq;
                    v.pp = tp; v.pq = tmpq;
                    // high/low-quality read split by goodq threshold (Java: tmpq >= goodq)
                    if (tmpq >= cfg_.goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
                    v.numberOfMismatches += nm;
                }
                qpos += len; rpe += len;
                break;
            }
            case 'D': {
                int p = rpos; // first deleted reference position
                if (p >= rlo && p <= rhi) {
                    std::string dels;
                    for (int i = 0; i < len; ++i) dels += ref_.at(rpos + i);
                    std::string sig = "-" + std::to_string(len);
                    out.positionToDeletionCount[p][sig]++;
                    int q = qpos < (int)bseq.size() ? bqual[qpos] : 30;
                    int tp = foldPos(rpe);
                    Variation& v = out.nonInsertionVariants[p][sig];
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += q;
                    v.meanMappingQuality += mapq;
                    v.numberOfMismatches += nm;
                    v.highQualityReadsCount++;
                }
                // addVariationForDeletion (CigarParser.java:1791): "increase coverage count for
                // reference bases missing from the read" -- a deletion read counts toward total
                // position coverage at every deleted base, so Depth at the deletion includes it.
                for (int i = 0; i < len; ++i)
                    if (rpos + i >= rlo && rpos + i <= rhi) out.refCoverage[rpos + i]++;
                rpos += len;
                break;
            }
            case 'N':
                rpos += len; // N (intron)
                break;
            case 'S': {
                // Faithful port of processSoftClip's mis-softclip re-matching + consensus storage
                // (chimeric SEED/SA detection is omitted; conf.chimeric defaults off and the seed
                // map is not built). The re-matching converts soft-clipped bases that actually match
                // the reference into reference-allele counts + coverage, which affects depth/AF near
                // clips; the remaining high-quality bases are stored as a soft-clip consensus for
                // realignment.
                bool isFivePrime = (k == 0);
                bool isThreePrime = (k == cigv.size() - 1);
                if (isFivePrime) {
                    int st = rpos;   // aligned start (VarDict 'start' == position)
                    int el = len;    // remaining soft-clip length
                    while (el - 1 >= 0 && st - 1 > 0 && ref_.has(st - 1) &&
                           bseq[el - 1] == ref_.at(st - 1) && bqual[el - 1] > 10) {
                        Variation& v = out.nonInsertionVariants[st - 1][std::string(1, ref_.at(st - 1))];
                        addCnt(v, reverse, el, bqual[el - 1], mapq, nm, cfg_.goodq);
                        if (st - 1 >= rlo && st - 1 <= rhi) out.refCoverage[st - 1]++;
                        st--; el--;
                    }
                    // Store high-quality remaining soft-clip as 5' consensus.
                    if (el > 0) {
                        int lowq = 0, hq = 0, sumq = 0;
                        for (int si = el - 1; si >= 0; --si) {
                            if (bseq[si] == 'N') break;
                            int bq = bqual[si];
                            if (bq <= 12) lowq++;
                            if (lowq > 1) break;
                            sumq += bq; hq++;
                        }
                        if (hq >= 1 && hq > lowq && st >= rlo && st <= rhi) {
                            Sclip& sc = out.softClips5End[st];
                            for (int si = el - 1; el - si <= hq; --si) {
                                char ch = bseq[si];
                                int idx = el - 1 - si;
                                sc.nt[idx][ch]++;
                                addCnt(getVariationFromSeq(sc, idx, ch),
                                       reverse, si - (el - hq), bqual[si], mapq, nm, cfg_.goodq);
                            }
                            addCnt(sc, reverse, len, (double)sumq / hq, mapq, nm, cfg_.goodq);
                        }
                    }
                } else if (isThreePrime) {
                    int st = rpos;    // current reference position
                    int qp = qpos;    // current query position
                    int el = len;
                    int rpeLocal = rpe;
                    while (qp < c.l_qseq && ref_.has(st) &&
                           bseq[qp] == ref_.at(st) && bqual[qp] > 10) {
                        Variation& v = out.nonInsertionVariants[st][std::string(1, ref_.at(st))];
                        addCnt(v, reverse, rlen - rpeLocal, bqual[qp], mapq, nm, cfg_.goodq);
                        if (st >= rlo && st <= rhi) out.refCoverage[st]++;
                        qp++; st++; el--; rpeLocal++;
                    }
                    if (c.l_qseq - qp > 0) {
                        int lowq = 0, hq = 0, sumq = 0;
                        for (int si = 0; si < el; ++si) {
                            if (bseq[qp + si] == 'N') break;
                            int bq = bqual[qp + si];
                            if (bq <= 12) lowq++;
                            if (lowq > 1) break;
                            sumq += bq; hq++;
                        }
                        if (hq >= 1 && hq > lowq && st >= rlo && st <= rhi) {
                            Sclip& sc = out.softClips3End[st];
                            for (int si = 0; si < hq; ++si) {
                                char ch = bseq[qp + si];
                                sc.nt[si][ch]++;
                                addCnt(getVariationFromSeq(sc, si, ch),
                                       reverse, hq - si, bqual[qp + si], mapq, nm, cfg_.goodq);
                            }
                            addCnt(sc, reverse, len, (double)sumq / hq, mapq, nm, cfg_.goodq);
                        }
                    }
                }
                qpos += len;
                break;
            }
            case 'H':
            case 'P':
            default:
                break;
            }
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(it);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp);
    return true;
}

} // namespace vardict
