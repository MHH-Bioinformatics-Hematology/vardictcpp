#pragma once
#include <string>
#include <vector>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// Emit the somatic (paired) 55-column rows for one region from the per-position candidate variants of
// the tumor (BAM1) and normal (BAM2) samples. Mirrors modes/SomaticMode +
// postprocessmodules/SomaticPostProcessModule (accept / callingForBothSamples / callingForOneSample).
void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& tumor,
                         const std::vector<SomaticPosition>& normal);

// Header for somatic mode (the two-sample 55-column layout); mirrors SomaticMode.printHeader.
void printSomaticHeader(std::FILE* out);

} // namespace vardict
