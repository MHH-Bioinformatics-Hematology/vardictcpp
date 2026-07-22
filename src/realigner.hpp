#pragma once
// Ports the parts of modules/VariationRealigner.java used by the simple pipeline.
#include "cigar_parser.hpp"
#include "reference.hpp"
#include "config.hpp"
#include "region.hpp"

namespace vardict {

// Merge adjacent SNVs / shorter MNVs into the full MNP they belong to (VariationRealigner.adjustMNP):
// counts move into the MNP, the partial variants are removed, and coverage is adjusted.
void adjustMNP(VariationData& vd, Reference& ref, const Config& cfg, const Region& region);

} // namespace vardict
