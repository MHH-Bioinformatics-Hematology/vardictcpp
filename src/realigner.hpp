#pragma once
// Ports the parts of modules/VariationRealigner.java used by the simple pipeline.
#include <functional>
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

// `reload(mstart,mend)` re-parses reads over a far breakpoint [mstart-200,mend+200] into vd (coverage
// only, no SV clusters), mirroring the partialPipeline reload that lets low-VAF SVs be AF-filtered.
using SVReloadFn = std::function<void(int, int)>;

// Realign large deletions inferred from soft-clip consensus breakpoints (findbp path). When the
// realigned breakpoint lands outside the current region (bp < region.start for a 5' clip, bp >
// region.end for a 3' clip), `reload` re-reads the coverage at the breakpoint so the deletion's
// AF reflects the true (usually high) depth there and is filtered exactly as in VarDict.
void realignlgdel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength, const SVReloadFn& reload);
void realignlgins30(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void adjSNV(VariationData& vd, Reference& ref);

// Structural variants (discordant-pair deletions). filterSVStructures collapses svfdel/svrdel mate
// clusters (VariationRealigner.filterSV) and must run before realignment; findDELdisc emits the <DEL>
// breakpoint variations (StructuralVariantsProcessor.findDELdisc) and runs after realignment.
void filterSVStructures(VariationData& vd, int maxReadLength);
void findDELdisc(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);

// StructuralVariantsProcessor.findINV: pair-assisted <INV> caller over the discordant same-orientation
// INV clusters. Runs after realignment, BEFORE findsv (so the split-read path skips folded soft clips).
void findINV(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
             int maxReadLength, const SVReloadFn& reload);
// StructuralVariantsProcessor.findsv: split-read SVs from soft clips. Only the candidate-inversion
// path is emitted here (<INV>); runs after realignment, before findDELdisc (findAllSVs order).
void findsv(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void realignlgins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);

} // namespace vardict
