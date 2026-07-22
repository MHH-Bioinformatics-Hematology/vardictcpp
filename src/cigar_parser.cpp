#include "cigar_parser.hpp"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <cstring>
#include <stdexcept>
#include <string>

namespace vardict {

static inline char baseChar(const bam1_t* b, int qpos) {
    static const char code[] = "=ACMGRSVTWYHKDBN";
    return code[bam_seqi(bam_get_seq(b), qpos)];
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

        bool reverse = (c.flag & BAM_FREVERSE) != 0;
        int mapq = c.qual;

        // NM tag (edit distance) if present.
        double nm = 0;
        if (uint8_t* aux = bam_aux_get(b, "NM")) nm = (double)bam_aux2i(aux);

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

        const uint32_t* cig = bam_get_cigar(b);
        const uint8_t* qual = bam_get_qual(b);
        int rpos = c.pos + 1;   // 1-based reference position of current op
        int qpos = 0;           // 0-based query offset (includes soft-clip)

        // VarDict read-position convention: position within the aligned read (M+I only, excluding
        // soft-clip), folded to the distance from the nearest read end. Precompute the aligned length.
        int rlen = 0;           // readLengthIncludeMatchingAndInsertions
        for (uint32_t k = 0; k < c.n_cigar; ++k) {
            int op = bam_cigar_op(cig[k]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CINS)
                rlen += bam_cigar_oplen(cig[k]);
        }
        int rpe = 0;            // readPositionExcludingSoftClipped
        auto foldPos = [&](int posExclSc) {
            return posExclSc < rlen - posExclSc ? posExclSc + 1 : rlen - posExclSc;
        };

        for (uint32_t k = 0; k < c.n_cigar; ++k) {
            int op = bam_cigar_op(cig[k]);
            int len = bam_cigar_oplen(cig[k]);
            switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for (int i = 0; i < len; ++i) {
                    int p = rpos + i;
                    int tp = foldPos(rpe + i);
                    if (p < rlo || p > rhi) continue;
                    char base = baseChar(b, qpos + i);
                    int q = qual[qpos + i];

                    out.refCoverage[p]++;
                    Variation& v = out.nonInsertionVariants[p][std::string(1, base)];
                    if (!v.pstd && v.pp != 0 && tp != v.pp) v.pstd = true;
                    if (!v.qstd && v.pq != 0 && (double)q != v.pq) v.qstd = true;
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += q;
                    v.meanMappingQuality += mapq;
                    v.numberOfMismatches += nm;
                    if (q >= cfg_.goodq) v.highQualityReadsCount++; else v.lowQualityReadsCount++;
                    v.pp = tp; v.pq = q;
                }
                rpos += len; qpos += len; rpe += len;
                break;
            case BAM_CINS: {
                int p = rpos - 1; // insertion anchored to preceding reference base (VarDict convention)
                if (p >= rlo && p <= rhi) {
                    std::string sig = "+";
                    double qsum = 0;
                    for (int i = 0; i < len; ++i) { sig += baseChar(b, qpos + i); qsum += qual[qpos + i]; }
                    int tp = foldPos(rpe);
                    Variation& v = out.insertionVariants[p][sig];
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += qsum / (len ? len : 1);
                    v.meanMappingQuality += mapq;
                    v.numberOfMismatches += nm;
                    v.highQualityReadsCount++;
                }
                qpos += len; rpe += len;
                break;
            }
            case BAM_CDEL: {
                int p = rpos; // first deleted reference position
                if (p >= rlo && p <= rhi) {
                    std::string dels;
                    for (int i = 0; i < len; ++i) dels += ref_.at(rpos + i);
                    std::string sig = "-" + std::to_string(len);
                    int q = qpos < c.l_qseq ? qual[qpos] : 30;
                    int tp = foldPos(rpe);
                    Variation& v = out.nonInsertionVariants[p][sig];
                    v.varsCount++;
                    v.incDir(reverse);
                    v.meanPosition += tp;
                    v.meanQuality += q;
                    v.meanMappingQuality += mapq;
                    v.numberOfMismatches += nm;
                    v.highQualityReadsCount++;
                    // deletion consumes reference coverage span for depth accounting
                    for (int i = 0; i < len; ++i) { int pp = rpos + i; if (pp >= rlo && pp <= rhi) out.refCoverage[pp]++; }
                }
                rpos += len;
                break;
            }
            case BAM_CREF_SKIP:
                rpos += len; // N (intron)
                break;
            case BAM_CSOFT_CLIP:
                qpos += len; // soft clip: consumes query only (consensus/realign not ported yet)
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
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
