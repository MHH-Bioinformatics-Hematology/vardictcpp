#pragma once
#include <map>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <htslib/sam.h>
#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "variation.hpp"

namespace vardict {

// Reusable BAM handle: opens the file, loads the index, and reads the header ONCE, then serves many
// per-region iterator queries. VarDictJava (htsjdk) keeps a single SamReader open for the whole run;
// the previous port re-ran sam_open/sam_index_load/sam_hdr_read for every region, so loading and
// parsing the (potentially large) .bai on each of ~20k regions dominated the wall time. One per worker
// thread — htslib file/index/header handles are not shared across threads.
class BamReader {
public:
    explicit BamReader(const std::string& bamPath) {
        fp_ = sam_open(bamPath.c_str(), "r");
        if (!fp_) throw std::runtime_error("cannot open BAM " + bamPath);
        idx_ = sam_index_load(fp_, bamPath.c_str());
        if (!idx_) throw std::runtime_error("cannot load BAM index for " + bamPath + " (run `samtools index`)");
        hdr_ = sam_hdr_read(fp_);
        if (!hdr_) throw std::runtime_error("cannot read BAM header");
    }
    ~BamReader() {
        if (hdr_) bam_hdr_destroy(hdr_);
        if (idx_) hts_idx_destroy(idx_);
        if (fp_) sam_close(fp_);
    }
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;
    samFile*    fp()  const { return fp_; }
    hts_idx_t*  idx() const { return idx_; }
    bam_hdr_t*  hdr() const { return hdr_; }
private:
    samFile*   fp_  = nullptr;
    hts_idx_t* idx_ = nullptr;
    bam_hdr_t* hdr_ = nullptr;
};

// Per-region variation data (mirrors data/scopedata/VariationData.java, counting subset).
struct VariationData {
    std::map<int, VarMap> nonInsertionVariants; // position -> allele -> Variation (ordered: drives emit)
    // These position-keyed maps are only ever accessed by key (never iterated in ascending order), or
    // are collected and std::sort'ed into a total order before use, so an unordered_map is byte-identical
    // and removes red-black-tree descents from the per-base counting hotpath.
    std::unordered_map<int, VarMap> insertionVariants;    // position -> "+SEQ" -> Variation
    std::unordered_map<int, int>    refCoverage;          // position -> total coverage
    std::unordered_map<int, std::map<std::string,int>> mnp; // position -> MNV description -> count
    std::unordered_map<int, std::map<std::string,int>> positionToInsertionCount; // pos -> "+SEQ" -> count
    std::unordered_map<int, std::map<std::string,int>> positionToDeletionCount;  // pos -> "-N" -> count
    std::map<int, Sclip> softClips5End;         // 5' soft-clip consensus per position
    std::map<int, Sclip> softClips3End;         // 3' soft-clip consensus per position
    // Structural-variant clusters (data/SVStructures.java): discordant read-pair deletion clusters,
    // forward (svfdel) / reverse (svrdel), plus their rolling right edges. Only the DEL discordant
    // path is collected/processed (findDELdisc); DUP/INV/fusion clusters are not built yet.
    std::vector<Sclip> svfdel, svrdel;
    int  svdelfend = 0, svdelrend = 0;
    std::map<int, SVInfo> svInfoAt;             // position -> SV marker (pairs/splits/clusters)
    int  maxReadLength = 0;
    int  chrLen = 0;            // length of the region's contig (for breakpoint bounds)
    long totalReads = 0;
    long dupReads = 0;
    double duprate() const { return totalReads > 0 ? (double)dupReads / (double)(totalReads + dupReads) : 0.0; }
};

// Reads the BAM over [region.start,region.end] (+halo) and accumulates per-position variation
// counts by walking each record's CIGAR. This is the memory-streaming heart of the port: records
// are consumed one at a time (htslib iterator), never held in a list. Ports CigarParser.java's
// per-base counting for M/=/X, I and D operations. Soft-clip consensus / realignment / SV are not
// yet implemented (documented in README) and so indel realignment adjustments are absent.
class CigarParser {
public:
    CigarParser(const Config& cfg, Reference& ref, BamReader& bam) : cfg_(cfg), ref_(ref), bam_(bam) {}

    // Fills `out` for the given region. Returns false on I/O error.
    bool process(const Region& region, VariationData& out);

private:
    const Config& cfg_;
    Reference& ref_;
    BamReader& bam_;
};

} // namespace vardict
