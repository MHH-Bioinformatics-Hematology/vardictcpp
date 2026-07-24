#pragma once
#include <string>
#include <vector>
#include <map>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// Per-amplicon, per-position variant set (mirrors variations/Vars.java as consumed by the amplicon
// post-processor). `variants` are the non-reference candidates (sorted as ToVarsBuilder.sortVariants,
// each carrying its good flag + description string + vartype). `hasRef`/`refTotalCov` stand in for the
// reference variant (only its coverage is consulted by AmpliconPostProcessModule).
struct AmpVars {
    std::vector<Variant> variants;
    bool hasRef = false;
    int  refTotalCov = 0;
};

// RegionBuilder.buildAmpRegions: parse an 8-column amplicon BED into segments. Regions are grouped by
// chromosome and merged into a segment while their insert intervals overlap/touch; a gap (or new
// chromosome) starts a new segment. AMP_BED_ROW_FORMAT = chr 0, start 1, end 2, gene 3, insertStart 6,
// insertEnd 7. numberNucleotideToExtend is NOT applied (unlike the plain BED path).
std::vector<std::vector<Region>> buildAmpRegions(const std::vector<std::string>& segRaws, const Config& cfg);

// Build the per-position AmpVars map for one amplicon region from its region-processing result.
std::map<int, AmpVars> buildAmpVars(const Config& cfg, const Region& region,
                                    const VariationData& vd, Reference& ref);

// AmpliconPostProcessModule.process for one segment: for every insert position covered by the segment's
// amplicons, pick the best call(s) across amplicons, flag amplicon bias, and append the 38-column rows.
void appendAmpliconSegment(std::string& out, const Config& cfg, const std::vector<Region>& regions,
                           const std::vector<std::map<int, AmpVars>>& vars);

// Amplicon-mode header (38 columns; AmpliconMode.printHeader).
void printAmpliconHeader(std::FILE* out);

} // namespace vardict
