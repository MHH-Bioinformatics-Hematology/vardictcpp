#pragma once
#include <map>
#include <string>
#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "variation.hpp"

namespace vardict {

// Per-region variation data (mirrors data/scopedata/VariationData.java, counting subset).
struct VariationData {
    std::map<int, VarMap> nonInsertionVariants; // position -> allele -> Variation
    std::map<int, VarMap> insertionVariants;    // position -> "+SEQ" -> Variation
    std::map<int, int>    refCoverage;          // position -> total coverage
    std::map<int, std::map<std::string,int>> mnp; // position -> MNV description -> count
    std::map<int, std::map<std::string,int>> positionToInsertionCount; // pos -> "+SEQ" -> count
    std::map<int, std::map<std::string,int>> positionToDeletionCount;  // pos -> "-N" -> count
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
    CigarParser(const Config& cfg, Reference& ref) : cfg_(cfg), ref_(ref) {}

    // Fills `out` for the given region. Returns false on I/O error.
    bool process(const Region& region, VariationData& out);

private:
    const Config& cfg_;
    Reference& ref_;
};

} // namespace vardict
