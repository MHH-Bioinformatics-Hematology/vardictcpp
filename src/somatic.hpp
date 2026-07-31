#pragma once
#include <string>
#include <vector>
#include <functional>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// combineAnalysis callback: runs the FULL pipeline over region2 = {chr, vstart-maxRL, vend+maxRL} on a
// MERGED bam1+bam2 read stream and looks up the variant at `position` whose descriptionString == desc.
// Returns true (and fills `out`) if that merged variant exists, false otherwise. Implemented in main.cpp
// (it captures both BamReaders, the Reference, runPipeline and callVariantsSomatic); consumed by
// SomaticPostProcessModule.combineAnalysis in somatic.cpp.
using CombineFn = std::function<bool(int vstart, int vend, int position,
                                     const std::string& desc, Variant& out)>;

// Emit the somatic (paired) 55-column rows for one region from the per-position candidate variants of
// the tumor (BAM1) and normal (BAM2) samples. Mirrors modes/SomaticMode +
// postprocessmodules/SomaticPostProcessModule (accept / callingForBothSamples / callingForOneSample).
// `combine` performs the merged-BAM refinement (combineAnalysis); pass an empty std::function to disable.
void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& tumor,
                         const std::vector<SomaticPosition>& normal,
                         const CombineFn& combine = CombineFn());

// Header for somatic mode (the two-sample 55-column layout); mirrors SomaticMode.printHeader.
void printSomaticHeader(std::FILE* out);

} // namespace vardict
