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

// Realign insertions/deletions: attribute nearby mismatch SNVs and soft-clip consensus to the
// indel (removing the spurious SNVs), and merge duplicate representations (VariationRealigner).
void realignins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void realigndel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);

// Realign large deletions inferred from soft-clip consensus breakpoints (findbp path).
void realignlgdel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void realignlgins30(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void adjSNV(VariationData& vd, Reference& ref);
void realignlgins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);

} // namespace vardict
