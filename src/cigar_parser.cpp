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
#include <vector>
#include <algorithm>

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

// CigarParser.isBEGIN_ATGC_AMP_ATGCs_END: true iff `s` is one A/T/G/C base, then '&', then one or
// more A/T/G/C bases (a pure MNV descriptor "X&YZ..."). Deletion ("-N&...") / insertion ("+...&...")
// descriptors do NOT match, so they are never recorded as MNVs.
static inline bool isMnpDesc(const std::string& s) {
    if (s.size() <= 2) return false;
    auto isATGC = [](char c) { return c == 'A' || c == 'T' || c == 'G' || c == 'C'; };
    if (s[1] != '&' || !isATGC(s[0])) return false;
    for (size_t i = 2; i < s.size(); ++i) if (!isATGC(s[i])) return false;
    return true;
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

// CigarParser.prepareSVStructuresForAnalysis (same-chromosome branch). Classifies a discordant read
// pair by orientation into a DEL / DUP / INV cluster and records its Mate, opening a new cluster when
// the read is > MINSVCDIST*maxReadLength past the last one of that kind. The cross-orientation
// adddisccnt() bumps are replicated exactly (Java 2056-2203) so every cluster's disc matches VarDict;
// DEL feeds findDELdisc, INV feeds findINV. DUP clusters are built only for that disc bookkeeping.
static void prepareSVStructures(const bam1_t* b, const Cig& cigv, int start,
                                const std::vector<int>& bqual, int lqseq, bool reverse,
                                double nm, VariationData& out, const Config& cfg) {
    const bam1_core_t& c = b->core;
    if (c.mtid != c.tid) return; // getMateReferenceName == "=" (same chr); inter-chr not ported
    const int maxRL = out.maxReadLength;
    const double CDIST = Config::MINSVCDIST * maxRL;
    const int MIN_D = 75;
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
    const double qAtBase = bqual[Config::MINMAPBASE];
    const double Qmean = c.qual;
    const double pmean = maxRL / 2.0;
    auto openIfNeeded = [](std::vector<Sclip>& lst, int startv, int fend, double cdist) {
        if (lst.empty() || startv - fend > cdist) { Sclip sc; sc.varsCount = 0; lst.push_back(sc); }
    };

    if (readDirNum * mateDirNum == -1 && (mlen * readDirNum) > 0 && lqseq > Config::MINMAPBASE) {
        // deletion candidate
        mlen = mateStart > start ? (long)mend - start : (long)end - mateStart;
        if (std::labs(mlen) > (long)cfg.INSSIZE + (long)cfg.INSSTDAMT * cfg.INSSTD) {
            if (readDirNum == 1) {
                openIfNeeded(out.svfdel, start, out.svdelfend, CDIST);
                addSVMate(out.svfdel.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft3, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svdelfend = end;
            } else {
                openIfNeeded(out.svrdel, start, out.svdelrend, CDIST);
                addSVMate(out.svrdel.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft5, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svdelrend = end;
            }
            if (!out.svfdel.empty() && std::abs(start - out.svdelfend) <= CDIST) out.svfdel.back().disc++;
            if (!out.svrdel.empty() && std::abs(start - out.svdelrend) <= CDIST) out.svrdel.back().disc++;
            if (!out.svfdup.empty() && std::abs(start - out.svdupfend) <= MIN_D) out.svfdup.back().disc++;
            if (!out.svrdup.empty() && std::abs(start - out.svduprend) <= MIN_D) out.svrdup.back().disc++;
            if (!out.svfinv5.empty() && std::abs(start - out.svinvfend5) <= MIN_D) out.svfinv5.back().disc++;
            if (!out.svrinv5.empty() && std::abs(start - out.svinvrend5) <= MIN_D) out.svrinv5.back().disc++;
            if (!out.svfinv3.empty() && std::abs(start - out.svinvfend3) <= MIN_D) out.svfinv3.back().disc++;
            if (!out.svrinv3.empty() && std::abs(start - out.svinvrend3) <= MIN_D) out.svrinv3.back().disc++;
        }
    } else if (readDirNum * mateDirNum == -1 && readDirNum * mlen < 0 && lqseq > Config::MINMAPBASE) {
        // duplication candidate
        if (readDirNum == 1) {
            openIfNeeded(out.svfdup, start, out.svdupfend, CDIST);
            addSVMate(out.svfdup.back(), start, end, mateStart, mend, readDirNum, totalLen,
                      (int)mlen, soft3, pmean, qAtBase, Qmean, nm, cfg.goodq);
            out.svdupfend = end;
        } else {
            openIfNeeded(out.svrdup, start, out.svduprend, CDIST);
            addSVMate(out.svrdup.back(), start, end, mateStart, mend, readDirNum, totalLen,
                      (int)mlen, soft5, pmean, qAtBase, Qmean, nm, cfg.goodq);
            out.svduprend = end;
        }
        if (!out.svfdup.empty() && std::abs(start - out.svdupfend) <= CDIST) out.svfdup.back().disc++;
        if (!out.svrdup.empty() && std::abs(start - out.svduprend) <= CDIST) out.svrdup.back().disc++;
        if (!out.svfdel.empty() && std::abs(start - out.svdelfend) <= MIN_D) out.svfdel.back().disc++;
        if (!out.svrdel.empty() && std::abs(start - out.svdelrend) <= MIN_D) out.svrdel.back().disc++;
        if (!out.svfinv5.empty() && std::abs(start - out.svinvfend5) <= MIN_D) out.svfinv5.back().disc++;
        if (!out.svrinv5.empty() && std::abs(start - out.svinvrend5) <= MIN_D) out.svrinv5.back().disc++;
        if (!out.svfinv3.empty() && std::abs(start - out.svinvfend3) <= MIN_D) out.svfinv3.back().disc++;
        if (!out.svrinv3.empty() && std::abs(start - out.svinvrend3) <= MIN_D) out.svrinv3.back().disc++;
    } else if (readDirNum * mateDirNum == 1 && lqseq > Config::MINMAPBASE) {
        // inversion candidate (read and mate same orientation)
        if (readDirNum == 1 && mlen != 0) {
            if (mlen < -3L * maxRL) {
                openIfNeeded(out.svfinv3, start, out.svinvfend3, CDIST);
                addSVMate(out.svfinv3.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft3, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svinvfend3 = end; out.svfinv3.back().disc++;
            } else if (mlen > 3L * maxRL) {
                openIfNeeded(out.svfinv5, start, out.svinvfend5, CDIST);
                addSVMate(out.svfinv5.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft3, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svinvfend5 = end; out.svfinv5.back().disc++;
            }
        } else if (mlen != 0) {
            if (mlen < -3L * maxRL) {
                openIfNeeded(out.svrinv3, start, out.svinvrend3, CDIST);
                addSVMate(out.svrinv3.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft5, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svinvrend3 = end; out.svrinv3.back().disc++;
            } else if (mlen > 3L * maxRL) {
                openIfNeeded(out.svrinv5, start, out.svinvrend5, CDIST);
                addSVMate(out.svrinv5.back(), start, end, mateStart, mend, readDirNum, totalLen,
                          (int)mlen, soft5, pmean, qAtBase, Qmean, nm, cfg.goodq);
                out.svinvrend5 = end; out.svrinv5.back().disc++;
            }
        }
        if (mlen != 0) {
            if (!out.svfdel.empty() && (start - out.svdelfend) <= MIN_D) out.svfdel.back().disc++;
            if (!out.svrdel.empty() && (start - out.svdelrend) <= MIN_D) out.svrdel.back().disc++;
            if (!out.svfdup.empty() && (start - out.svdupfend) <= MIN_D) out.svfdup.back().disc++;
            if (!out.svrdup.empty() && (start - out.svduprend) <= MIN_D) out.svrdup.back().disc++;
        }
    }
}

bool CigarParser::process(const Region& region, VariationData& out, bool reloadMode) {
    // File/index/header are opened once per worker thread (BamReader) and reused across all regions.
    samFile*    fp  = bam_.fp();
    hts_idx_t*  idx = bam_.idx();
    bam_hdr_t*  hdr = bam_.hdr();

    int tid = bam_name2id(hdr, region.chr.c_str());
    if (tid < 0) { // try toggling chr prefix
        std::string alt = (region.chr.rfind("chr", 0) == 0) ? region.chr.substr(3) : "chr" + region.chr;
        tid = bam_name2id(hdr, alt.c_str());
    }
    if (tid < 0) { return true; }

    out.chrLen = hdr->target_len[tid];
    // Fetch reads overlapping the region (htslib is 0-based, end-exclusive).
    hts_itr_t* it = sam_itr_queryi(idx, tid, region.start - 1, region.end);
    bam1_t* b = bam_init1();

    const int rlo = region.start, rhi = region.end;
    // Simple duplicate accounting for the duprate column (start+cigar signature per alignment start).
    long lastStart = -1;
    std::map<std::string, int> dupKeys;

    // Per-read scratch buffers, reused across reads: reads are near-uniform length, so after the first
    // record the string/vectors keep their capacity and no longer heap-allocate each iteration.
    std::string bseq;
    std::vector<int> bqual;
    Cig cigv;

    while (sam_itr_next(fp, it, b) >= 0) {
        const bam1_core_t& c = b->core;
        if (c.flag & BAM_FUNMAP) continue;
        if (cfg_.samFilterFlag != 0 && (c.flag & cfg_.samFilterFlag)) continue;
        if ((double)c.qual < cfg_.mapqMin) continue;
        if (c.n_cigar == 0 || c.l_qseq == 0) continue;
        // Ignore supplementary alignments (they would skew coverage), as VarDict does.
        if (c.flag & BAM_FSUPPLEMENTARY) continue;
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
        bseq.assign(c.l_qseq, 'N');
        bqual.resize(c.l_qseq);
        for (int j = 0; j < c.l_qseq; ++j) { bseq[j] = CODE[bam_seqi(rawseq, j)]; bqual[j] = rawqual[j]; }
        const uint32_t* rawcig = bam_get_cigar(b);
        static const char OPS[] = "MIDNSHP=X";
        cigv.clear();
        for (uint32_t k = 0; k < c.n_cigar; ++k) cigv.push_back({(int)bam_cigar_oplen(rawcig[k]), OPS[bam_cigar_op(rawcig[k])]});
        int rpos = c.pos + 1;   // 1-based reference position of current op
        if (cfg_.performLocalRealignment) modifyCigar(rpos, cigv, bseq, bqual, ref_, out.maxReadLength, cfg_);
        // Ignore reads soft-clipped at both ends where the leading clip is 10-99 bp and the trailing
        // clip is >= 10 bp (VarDict pattern ^\d\dS.*\d\dS$: two-digit leading S, >=2-digit trailing S).
        // VarDict (CigarParser) applies this to the POST-modifyCigar CIGAR, so a read whose leading
        // insertion/short-match is turned into a soft-clip by modifyCigar can newly qualify.
        if (!cigv.empty()) {
            int lead = (cigv.front().second == 'S') ? cigv.front().first : 0;
            int tail = (cigv.back().second == 'S') ? cigv.back().first : 0;
            if (lead >= 10 && lead <= 99 && tail >= 10) continue;
        }
        // maxReadLength is updated AFTER modifyCigar (CigarParser l.307-309: getSoftClippedLength of
        // the MODIFIED cigar), so the CigarModifier chimeric-clip check for THIS read sees only the
        // PREVIOUS reads' max -- a read cannot use its own length to justify dropping its own soft
        // clip via the `abs(pos - seed) < 2*maxReadLength` gate. Value = M+I+S of the modified cigar.
        {
            int solen = 0;
            for (auto& e : cigv)
                if (e.second=='M'||e.second=='I'||e.second=='S'||e.second=='='||e.second=='X') solen += e.first;
            if (solen > out.maxReadLength) out.maxReadLength = solen;
        }
        const int readStart = rpos;  // CigarModifier-adjusted read alignment start (Java's `position`)

        // Structural-variant discordant-pair collection (CigarParser dispatch at 323-329): skip
        // paired reads whose mate is unmapped (potential insertion, not ported); otherwise, for
        // MAPQ>10 reads, record possible SV deletion clusters from the (post-modifyCigar) alignment.
        // In reloadMode (StructuralVariantsProcessor's partialPipeline) the reload's CigarParser is
        // constructed with a throwaway SVStructures, so its SV-cluster bumps are discarded; we simply
        // skip them here. Coverage/variations/soft-clips still accumulate into the shared `out`.
        if (!cfg_.disableSV && !reloadMode) {
            bool paired = (c.flag & BAM_FPAIRED) != 0;
            bool mateUnmapped = (c.flag & BAM_FMUNMAP) != 0;
            if (paired && mateUnmapped) { /* potential insertion: not ported */ }
            else if (c.qual > 10) prepareSVStructures(b, cigv, rpos, bqual, c.l_qseq, reverse, nm, out, cfg_);
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

        // CigarParser 'offset': number of bases of the FOLLOWING matched segment already consumed by
        // a complex variant on the current segment (findOffset). Persists into the next M segment
        // (its inner loop starts at i=offset) and is reset to 0 on I/D/H/S segments. Reset per read.
        int carryOffset = 0;

        // CigarParser.findOffset: scan up to `cigarLen` bases of the next matched segment, extending an
        // offset while bases stay within conf.vext of a mismatch. Returns how many bases to fold into
        // the complex variant, the base sequence, per-base quality values, and the mismatch count; also
        // bumps refCoverage for the folded bases (as VarDict does).
        struct OffsetRes { int offset; std::string ss; std::vector<int> quals; int tnm; };
        auto findOffset = [&](int refPos, int readPos, int cigarLen) -> OffsetRes {
            OffsetRes r{0, "", {}, 0};
            int vsn = 0;
            for (int vi = 0; vsn <= cfg_.vext && vi < cigarLen; vi++) {
                if (readPos + vi >= (int)bseq.size()) break;
                if (bseq[readPos + vi] == 'N') break;
                if (bqual[readPos + vi] < cfg_.goodq) break;
                if (ref_.has(refPos + vi)) {
                    if (bseq[readPos + vi] != ref_.at(refPos + vi)) { r.offset = vi + 1; r.tnm++; vsn = 0; }
                    else vsn++;
                }
            }
            if (r.offset > 0) {
                r.ss = bseq.substr(readPos, r.offset);
                for (int osi = 0; osi < r.offset; osi++) {
                    r.quals.push_back(bqual[readPos + osi]);
                    int cp = refPos + osi;
                    if (cp >= rlo && cp <= rhi) out.refCoverage[cp]++;
                }
            }
            return r;
        };

        for (uint32_t k = 0; k < cigv.size(); ++k) {
            char op = cigv[k].second;
            int len = cigv[k].first;
            switch (op) {
            case 'M':
            case '=':
            case 'X': {
                // Faithful port of CigarParser's matching-part loop (per-base, mutating rpos/qpos/rpe
                // like VarDictJava's start/readPositionIncludingSoftClipped/readPositionExcludingSoftClipped).
                // Grows adjacent mismatches into MNVs ("A&CG"), and — the piece previously deferred —
                // bridges the trailing base(s) of this M-segment into a following Deletion ("-N&...",
                // with '^' two-insertions-ahead and a findOffset '&' tail) or Insertion ("+X&YZ"),
                // producing one complex variant instead of a split SNV + indel.
                int nmoff = 0;      // consecutive mismatches merged into MNVs (not double-counted as "other")
                int moffset = 0;    // findOffset carry into the next M segment
                int i = carryOffset;
                while (i < len) {
                    char ch1 = bseq[qpos];
                    if (ch1 == 'N') { rpos++; qpos++; rpe++; i++; continue; }
                    double q = bqual[qpos];
                    int qbases = 1;     // bases counted for reference coverage / position averaging
                    int qibases = 0;    // inserted bases counted only for quality averaging
                    std::string s(1, ch1);
                    std::string ss;
                    bool startWithDeletion = false;
                    int ddlen = 0;
                    // Grow MNV while the current base mismatches the reference and quality is good.
                    while ((rpos + 1) >= rlo && (rpos + 1) <= rhi && (i + 1) < len &&
                           q >= cfg_.goodq &&
                           ref_.has(rpos) && ref_.at(rpos) != bseq[qpos] && ref_.at(rpos) != 'N') {
                        if (bqual[qpos + 1] < cfg_.goodq + 5) break;
                        char nuc = bseq[qpos + 1];
                        if (nuc == 'N') break;
                        if (ref_.has(rpos + 1) && ref_.at(rpos + 1) == 'N') break;
                        if (!(ref_.has(rpos + 1) && ref_.at(rpos + 1) == nuc)) {  // next base also mismatches
                            ss += nuc; q += bqual[qpos + 1]; qbases++;
                            qpos++; rpe++; i++; rpos++;
                            nmoff++;
                        } else {                                  // bridge matching bases to next mismatch within vext
                            int ssn = 0;
                            for (int ssi = 1; ssi <= cfg_.vext; ssi++) {
                                if (i + 1 + ssi >= len) break;
                                if (qpos + 1 + ssi < (int)bseq.size() && ref_.has(rpos + 1 + ssi) &&
                                    bseq[qpos + 1 + ssi] != ref_.at(rpos + 1 + ssi)) { ssn = ssi + 1; break; }
                            }
                            if (ssn == 0) break;
                            if (bqual[qpos + ssn] < cfg_.goodq + 5) break;
                            for (int ssi = 1; ssi <= ssn; ssi++) { ss += bseq[qpos + ssi]; q += bqual[qpos + ssi]; qbases++; }
                            qpos += ssn; rpe += ssn; i += ssn; rpos += ssn;
                        }
                    }
                    if (!ss.empty()) s += "&" + ss;

                    // isCloserThenVextAndGoodBase: near the end of this M segment (within vext), with a
                    // mismatch (or grown MNV) at good quality and the next CIGAR op being cop.
                    auto closerGood = [&](char cop) -> bool {
                        if (k + 2 < cigv.size() && cigv[k + 2].second == 'H') return false;
                        return cfg_.performLocalRealignment && (len - i) <= cfg_.vext &&
                               (k + 1) < cigv.size() && cigv[k + 1].second == cop &&
                               ref_.has(rpos) &&
                               (!ss.empty() || bseq[qpos] != ref_.at(rpos)) &&
                               bqual[qpos] >= cfg_.goodq;
                    };

                    if (closerGood('D')) {
                        while (i + 1 < len) {                     // fold remaining M bases into s
                            s += bseq[qpos + 1]; q += bqual[qpos + 1]; qbases++;
                            i++; qpos++; rpe++; rpos++;
                        }
                        auto amp = s.find('&'); if (amp != std::string::npos) s.erase(amp, 1);
                        ddlen = cigv[k + 1].first;
                        s = "-" + std::to_string(ddlen) + "&" + s;
                        startWithDeletion = true;
                        k += 1;                                   // consume the D segment
                        if (k + 1 < cigv.size() && cigv[k + 1].second == 'I') {   // insertion two ahead
                            int n2 = cigv[k + 1].first;
                            s += "^" + bseq.substr(qpos + 1, n2);
                            for (int qi = 1; qi <= n2; qi++) { q += bqual[qpos + 1 + qi]; qibases++; }
                            qpos += n2; rpe += n2;
                            k += 1;                               // consume the I segment
                        }
                        if (k + 1 < cigv.size() && cigv[k + 1].second == 'M') {   // extend into next M
                            OffsetRes tpl = findOffset(rpos + ddlen + 1, qpos + 1, cigv[k + 1].first);
                            if (tpl.offset != 0) {
                                moffset = tpl.offset;
                                nmoff += tpl.tnm;
                                s += "&" + tpl.ss;
                                for (int qv : tpl.quals) { q += qv; qibases++; }
                            }
                        }
                    } else if (closerGood('I')) {
                        while (i + 1 < len) {                     // fold remaining M bases into s
                            s += bseq[qpos + 1]; q += bqual[qpos + 1]; qbases++;
                            i++; qpos++; rpe++; rpos++;
                        }
                        auto amp = s.find('&'); if (amp != std::string::npos) s.erase(amp, 1);
                        int n2 = cigv[k + 1].first;
                        s += bseq.substr(qpos + 1, n2);
                        s = s.substr(0, n2) + "&" + s.substr(n2);
                        s = "+" + s;
                        for (int qi = 1; qi <= n2; qi++) { q += bqual[qpos + 1 + qi]; qibases++; }
                        qpos += n2; rpe += n2;
                        k += 1;                                   // consume the I segment
                        qibases--; qbases++;                      // set the correct insertion anchor position
                    }

                    int pos = rpos - qbases + 1;                  // leftmost covered position
                    if (pos >= rlo && pos <= rhi && s.find('N') == std::string::npos) {
                        double qavg = q / (qbases + qibases);
                        int tp = rpe < rlen - rpe ? rpe + 1 : rlen - rpe;
                        bool isIns = (!s.empty() && s[0] == '+');
                        if (isIns) out.positionToInsertionCount[pos][s]++;
                        else if (isMnpDesc(s)) out.mnp[pos][s]++;
                        Variation& v = isIns ? out.insertionVariants[pos][s] : out.nonInsertionVariants[pos][s];
                        if (!v.pstd && v.pp != 0 && tp != v.pp) v.pstd = true;
                        if (!v.qstd && v.pq != 0 && qavg != v.pq) v.qstd = true;
                        v.varsCount++;
                        v.incDir(reverse);
                        v.meanPosition += tp;
                        v.meanQuality += qavg;
                        v.meanMappingQuality += mapq;
                        v.numberOfMismatches += nm - nmoff;
                        v.pp = tp; v.pq = qavg;
                        if (qavg >= cfg_.goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
                        // reference coverage for every base covered by this variation
                        int shift = (isIns && s.find('&') != std::string::npos) ? 1 : 0;
                        for (int qi = 1; qi <= qbases - shift; ++qi) {
                            int cp = rpos - qi + 1;
                            if (cp >= rlo && cp <= rhi) out.refCoverage[cp]++;
                        }
                        if (startWithDeletion) {
                            out.positionToDeletionCount[pos][s]++;
                            for (int qi = 1; qi < ddlen; qi++) {
                                int cp = rpos + qi;
                                if (cp >= rlo && cp <= rhi) out.refCoverage[cp]++;
                            }
                        }
                    }
                    if (startWithDeletion) rpos += ddlen;
                    rpos++; qpos++; rpe++;                        // advance past the current base (M: ref+read)
                    i++;
                }
                carryOffset = 0;
                if (moffset != 0) { carryOffset = moffset; qpos += moffset; rpos += moffset; rpe += moffset; }
                break;
            }
            case 'I': {
                carryOffset = 0;  // CigarParser resets offset before processing an insertion
                int p = rpos - 1; // insertion anchored to preceding reference base (VarDict convention)
                std::string ins;
                double qsum = 0;
                for (int i = 0; i < len; ++i) { ins += bseq[qpos + i]; qsum += bqual[qpos + i]; }

                // processInsertion's two complex-tail paths:
                //  - isInsertionOrDeletionWithNextMatched (I + short M(<=vext) + I/D + non-indel): fold
                //    the trailing M and following indel into the insertion via appendSegments,
                //    producing "+ins#matchedM^indel" (+ a findOffset "&ss" tail), then advance ref/read
                //    by the extra segments (multoffs/multoffp).
                //  - isNextMatched (I + M): fold leading mismatch(es) of the next M into "+ins&mismatch".
                // qcount tracks Java's qualityString length so the mean base quality (tmpq) matches.
                int offset = 0, nmoff = 0, multoffs = 0, multoffp = 0, qcount = len;
                std::string ssx;
                std::string desc = ins;               // descStringOfInsertionSegment
                bool insDelNext = cfg_.performLocalRealignment && (k + 2) < cigv.size()
                    && cigv[k + 1].first <= cfg_.vext && cigv[k + 1].second == 'M'
                    && (cigv[k + 2].second == 'I' || cigv[k + 2].second == 'D')
                    && ((k + 3) >= cigv.size() || (cigv[k + 3].second != 'I' && cigv[k + 3].second != 'D'));
                if (insDelNext) {
                    int mLen = cigv[k + 1].first, indelLen = cigv[k + 2].first, begin = qpos + len;
                    // appendSegments (isInsertion=true): "#" + matched M, then "^" + (I seq | D len).
                    desc += "#" + bseq.substr(begin, mLen);
                    for (int i = 0; i < mLen; i++) { qsum += bqual[begin + i]; qcount++; }
                    if (cigv[k + 2].second == 'I') {
                        desc += "^" + bseq.substr(begin + mLen, indelLen);
                        for (int i = 0; i < indelLen; i++) { qsum += bqual[begin + mLen + i]; qcount++; }
                    } else {
                        desc += "^" + std::to_string(indelLen);
                        qsum += bqual[begin + mLen]; qcount++;   // isInsertion=true: one quality char for D
                    }
                    multoffs += mLen + (cigv[k + 2].second == 'D' ? indelLen : 0);
                    multoffp += mLen + (cigv[k + 2].second == 'I' ? indelLen : 0);
                    if (k + 3 < cigv.size() && cigv[k + 3].second == 'M') {
                        OffsetRes tpl = findOffset(rpos + multoffs, qpos + len + multoffp, cigv[k + 3].first);
                        offset = tpl.offset; ssx = tpl.ss;
                        for (int qv : tpl.quals) { qsum += qv; qcount++; }   // nmoff NOT bumped here (Java)
                    }
                    k += 2;                             // consume the M and indel segments
                } else if (cfg_.performLocalRealignment && (k + 1) < cigv.size() && cigv[k + 1].second == 'M') {
                    int mlen = cigv[k + 1].first;
                    int mbeg = qpos + len;   // readPositionIncludingSoftClipped + cigarElementLength
                    int vsn = 0;
                    for (int vi = 0; vsn <= cfg_.vext && vi < mlen; vi++) {
                        if (mbeg + vi >= (int)bseq.size()) break;
                        if (bseq[mbeg + vi] == 'N') break;
                        if (bqual[mbeg + vi] < cfg_.goodq) break;
                        if (ref_.has(rpos + vi)) {
                            if (bseq[mbeg + vi] != ref_.at(rpos + vi)) { offset = vi + 1; nmoff++; vsn = 0; }
                            else vsn++;
                        }
                    }
                    if (offset != 0) {
                        ssx = bseq.substr(mbeg, offset);
                        for (int osi = 0; osi < offset; osi++) {
                            qsum += bqual[mbeg + osi]; qcount++;
                            int cp = rpos + osi;
                            if (cp >= rlo && cp <= rhi) out.refCoverage[cp]++;
                        }
                    }
                }

                if (offset > 0) desc += "&" + ssx;
                // adjInsPos only for pure-ATGC insertion descriptors (BEGIN_ATGC_END); complex
                // "#/^/&" descriptors keep the anchor at rpos-1.
                if (desc.find_first_not_of("ACGT") == std::string::npos && !desc.empty())
                    adjInsPos(p, desc, ref_);          // left-normalize insertion anchor in repeats
                std::string sig = "+" + desc;
                if (p >= rlo && p <= rhi && desc.find('N') == std::string::npos) {
                    int tp = foldPos(rpe);
                    // mean quality over the inserted segment plus any folded/appended bases (Java tmpq)
                    double tmpq = qsum / (qcount ? qcount : 1);
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
                    v.numberOfMismatches += nm - nmoff;

                    // Finding 3 (processInsertion subCnt): the anchor base (read[qpos-1], ref position p)
                    // was counted as a reference observation by the preceding M run; this read actually
                    // supports the insertion, so remove that contribution from the anchor ref allele when
                    // the anchor base matches the reference. refCoverage is left intact. The gate is the
                    // read's own alignment start (Java processInsertion: `insertionPosition > position`),
                    // NOT the region start -- otherwise an insertion at the first position of a region
                    // wrongly keeps the anchor ref count (leaving a residual RefFwd/RefRev read).
                    if (p > readStart && ref_.has(p) && bseq[qpos - 1] == ref_.at(p)) {
                        auto pit = out.nonInsertionVariants.find(p);
                        if (pit != out.nonInsertionVariants.end()) {
                            auto vit = pit->second.find(std::string(1, bseq[qpos - 1]));
                            if (vit != pit->second.end()) {
                                Variation& tv = vit->second;
                                double bq = bqual[qpos - 1];
                                tv.varsCount--;
                                tv.decDir(reverse);
                                tv.meanPosition -= tp;
                                tv.meanQuality -= bq;
                                tv.meanMappingQuality -= mapq;
                                tv.numberOfMismatches -= (nm - nmoff);
                                if (bq >= cfg_.goodq) tv.highQualityReadsCount--; else tv.lowQualityReadsCount--;
                            }
                        }
                    }
                    // Finding 4: insertion at the read edge (2nd CIGAR op, first op a soft/hard clip) --
                    // add one anchor reference observation + coverage so the insertion AF can't exceed 1.
                    if (k == 1 && !cigv.empty() && (cigv[0].second == 'S' || cigv[0].second == 'H') && ref_.has(p)) {
                        Variation& tt = out.nonInsertionVariants[p][std::string(1, ref_.at(p))];
                        tt.incDir(reverse);
                        tt.varsCount++;
                        tt.pstd = v.pstd; tt.qstd = v.qstd;
                        tt.meanPosition += tp;
                        tt.meanQuality += tmpq;
                        tt.meanMappingQuality += mapq;
                        tt.pp = tp; tt.pq = tmpq;
                        tt.numberOfMismatches += nm - nmoff;
                        out.refCoverage[p]++;
                    }
                }
                qpos += len + offset + multoffp; rpe += len + offset + multoffp; rpos += offset + multoffs;
                carryOffset = offset;   // next M segment starts past the folded bases
                break;
            }
            case 'D': {
                // Faithful port of CigarParser.processDeletion: the deletion description grows a
                // trailing complex tail when the following CIGAR op is matched/inserted, producing
                // "-N&ss" (D+M mismatch), "-N^ins" (D+I) or "-N#seg^..." (D+M+indel) instead of a
                // bare "-N" that would split into a separate SNV/indel downstream.
                carryOffset = 0;  // CigarParser resets offset before processing a deletion
                // skipIndelNextToIntron: deletions adjacent to an intron (N) are ignored (RNA-seq).
                if ((k + 1 < cigv.size() && cigv[k + 1].second == 'N') ||
                    (k > 0 && cigv[k - 1].second == 'N')) {
                    rpe += len;
                    break;
                }
                std::string descStr = "-" + std::to_string(len);  // descStringOfDeletedElement
                std::string ssAppend;                              // sequenceToAppendIfNextSegmentMatched
                std::vector<int> qualSeg;                          // qualityOfSegment (decoded values)
                int q1 = bqual[qpos - 1];                          // qualityOfLastSegmentBeforeDel
                int multoffs = 0, multoffp = 0, nmoff = 0, offset = 0;

                bool lr = cfg_.performLocalRealignment;
                // isInsertionOrDeletionWithNextMatched: D, short M(<=vext), then I/D, then non-indel.
                bool branchA = lr && k + 3 < cigv.size() &&
                               cigv[k + 1].first <= cfg_.vext && cigv[k + 1].second == 'M' &&
                               (cigv[k + 2].second == 'I' || cigv[k + 2].second == 'D') &&
                               cigv[k + 3].second != 'I' && cigv[k + 3].second != 'D';
                bool branchB = lr && k + 1 < cigv.size() && cigv[k + 1].second == 'I';
                bool branchC = lr && k + 1 < cigv.size() && cigv[k + 1].second == 'M';

                if (branchA) {
                    int mLen = cigv[k + 1].first, indelLen = cigv[k + 2].first, begin = qpos;
                    // appendSegments (isInsertion=false)
                    descStr += "#" + bseq.substr(begin, mLen);
                    for (int i = 0; i < mLen; i++) qualSeg.push_back(bqual[begin + i]);
                    if (cigv[k + 2].second == 'I') {
                        descStr += "^" + bseq.substr(begin + mLen, indelLen);
                        for (int i = 0; i < indelLen; i++) qualSeg.push_back(bqual[begin + mLen + i]);
                    } else {
                        descStr += "^" + std::to_string(indelLen);  // D two-ahead: no quality appended
                    }
                    multoffs += mLen + (cigv[k + 2].second == 'D' ? indelLen : 0);
                    multoffp += mLen + (cigv[k + 2].second == 'I' ? indelLen : 0);
                    if (k + 3 < cigv.size() && cigv[k + 3].second == 'M') {
                        int vsn = 0, tn = qpos + multoffp, ts = rpos + multoffs + len, seglen = cigv[k + 3].first;
                        for (int vi = 0; vsn <= cfg_.vext && vi < seglen; vi++) {
                            if (tn + vi >= (int)bseq.size() || bseq[tn + vi] == 'N') break;
                            if (bqual[tn + vi] < cfg_.goodq) break;
                            if (ref_.has(ts + vi) && ref_.at(ts + vi) == 'N') break;
                            if (ref_.has(ts + vi)) {
                                if (bseq[tn + vi] != ref_.at(ts + vi)) { offset = vi + 1; nmoff++; vsn = 0; }
                                else vsn++;
                            }
                        }
                        if (offset != 0) {
                            ssAppend += bseq.substr(tn, offset);
                            for (int i = 0; i < offset; i++) qualSeg.push_back(bqual[tn + i]);
                        }
                    }
                    k += 2;
                } else if (branchB) {
                    int insLen = cigv[k + 1].first;
                    descStr += "^" + bseq.substr(qpos, insLen);
                    for (int i = 0; i < insLen; i++) qualSeg.push_back(bqual[qpos + i]);
                    multoffp += insLen;
                    if (k + 2 < cigv.size() && cigv[k + 2].second == 'M') {
                        int mLen = cigv[k + 2].first, vsn = 0, tn = qpos + multoffp, ts = rpos + len;
                        for (int vi = 0; vsn <= cfg_.vext && vi < mLen; vi++) {
                            if (tn + vi >= (int)bseq.size() || bseq[tn + vi] == 'N') break;
                            if (bqual[tn + vi] < cfg_.goodq) break;
                            if (ref_.has(ts + vi)) {
                                if (ref_.at(ts + vi) == 'N') break;
                                if (bseq[tn + vi] != ref_.at(ts + vi)) { offset = vi + 1; nmoff++; vsn = 0; }
                                else vsn++;
                            }
                        }
                        if (offset != 0) {
                            ssAppend += bseq.substr(tn, offset);
                            for (int i = 0; i < offset; i++) qualSeg.push_back(bqual[tn + i]);
                        }
                    }
                    k += 1;
                } else if (branchC) {
                    int mLen = cigv[k + 1].first, vsn = 0;
                    for (int vi = 0; vsn <= cfg_.vext && vi < mLen; vi++) {
                        if (qpos + vi >= (int)bseq.size() || bseq[qpos + vi] == 'N') break;
                        if (bqual[qpos + vi] < cfg_.goodq) break;
                        if (ref_.has(rpos + len + vi)) {
                            if (ref_.at(rpos + len + vi) == 'N') break;
                            if (bseq[qpos + vi] != ref_.at(rpos + len + vi)) { offset = vi + 1; nmoff++; vsn = 0; }
                            else vsn++;
                        }
                    }
                    if (offset != 0) {
                        ssAppend += bseq.substr(qpos, offset);
                        for (int i = 0; i < offset; i++) qualSeg.push_back(bqual[qpos + i]);
                    }
                }

                if (offset > 0) descStr += "&" + ssAppend;
                // quality of first matched base after the deletion: best of q1 and q2.
                if (qpos + offset >= (int)bseq.size()) qualSeg.push_back(q1);
                else { int q2 = bqual[qpos + offset]; qualSeg.push_back(q1 > q2 ? q1 : q2); }

                // addVariationForDeletion (only when the deletion start is inside the region).
                if (rpos >= rlo && rpos <= rhi) {
                    out.positionToDeletionCount[rpos][descStr]++;
                    Variation& v = out.nonInsertionVariants[rpos][descStr];
                    int tp = foldPos(rpe);
                    double tmpq = 0; for (int qv : qualSeg) tmpq += qv; tmpq /= qualSeg.size();
                    if (!v.pstd && v.pp != 0 && tp != v.pp) v.pstd = true;
                    if (!v.qstd && v.pq != 0 && tmpq != v.pq) v.qstd = true;
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += tmpq;
                    v.meanMappingQuality += mapq;
                    v.pp = tp; v.pq = tmpq;
                    v.numberOfMismatches += nm - nmoff;
                    if (tmpq >= cfg_.goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
                    // increase coverage for reference bases missing from the read
                    for (int i = 0; i < len; ++i) out.refCoverage[rpos + i]++;
                }
                rpos += len + offset + multoffs;
                qpos += offset + multoffp;
                rpe += offset + multoffp;
                carryOffset = offset;
                break;
            }
            case 'N':
                // Intron: record the splice junction (CigarParser.processNotMatched). isGoodVar rejects
                // a Deletion whose coordinates equal a splice junction (it is an intron, not a deletion).
                out.splice.insert(std::to_string(rpos - 1) + "-" + std::to_string(rpos + len - 1));
                rpos += len;
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
                            // Java sclip5HighQualityProcessing passes the *remaining* soft-clip length
                            // (cigarElementLength, decremented by the mis-softclip match-extension loop
                            // above), not the original clip length, into the whole-Sclip addCnt.
                            addCnt(sc, reverse, el, (double)sumq / hq, mapq, nm, cfg_.goodq);
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
                            // Java sclip3HighQualityProcessing passes the *remaining* soft-clip length
                            // (cigarElementLength, decremented by the mis-softclip match-extension loop
                            // above), not the original clip length, into the whole-Sclip addCnt.
                            addCnt(sc, reverse, el, (double)sumq / hq, mapq, nm, cfg_.goodq);
                        }
                    }
                }
                qpos += len;
                carryOffset = 0;  // processSoftClip resets offset
                break;
            }
            case 'H':
                carryOffset = 0;  // hard-clip resets offset
                break;
            case 'P':
            default:
                break;
            }
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(it);
    return true;
}

} // namespace vardict
