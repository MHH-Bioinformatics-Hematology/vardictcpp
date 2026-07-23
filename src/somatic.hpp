#pragma once
#include <string>
#include <vector>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// Emit the somatic (paired) 55-column rows for one region. `positions` holds the per-position candidate
// variants for the (tumor, normal) pair; when the two BAMs are identical the same set is used for both
// samples. Mirrors modes/SomaticMode + postprocessmodules/SomaticPostProcessModule.
void appendSomaticRegion(std::string& out, const Config& cfg, const Region& region,
                         const std::vector<SomaticPosition>& positions);

// Header for somatic mode (57 fields incl. the two-sample layout); mirrors SomaticMode.printHeader.
void printSomaticHeader(std::FILE* out);

} // namespace vardict
