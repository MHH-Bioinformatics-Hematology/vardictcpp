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
// `bams` are the open BAM readers backing this pipeline (single-sample: one; merged combineAnalysis:
// both). realigndel needs them for the noPassingReads microhomology check (VariationRealigner l.617);
// pass empty to skip it (behaves as bams==null in Java).
void realigndel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength,
                const std::vector<BamReader*>& bams);

// `reload(mstart,mend)` re-parses reads over a far breakpoint [mstart-200,mend+200] into vd (coverage
// only, no SV clusters), mirroring the partialPipeline reload that lets low-VAF SVs be AF-filtered.
using SVReloadFn = std::function<void(int, int)>;

// Realign large deletions inferred from soft-clip consensus breakpoints (findbp path). When the
// realigned breakpoint lands outside the current region (bp < region.start for a 5' clip, bp >
// region.end for a 3' clip), `reload` re-reads the coverage at the breakpoint so the deletion's
// AF reflects the true (usually high) depth there and is filtered exactly as in VarDict.
void realignlgdel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength,
                  const SVReloadFn& reload, const std::vector<BamReader*>& bams);
void realignlgins30(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
void adjSNV(VariationData& vd, Reference& ref);

// Structural variants (discordant-pair deletions). filterSVStructures collapses svfdel/svrdel mate
// clusters (VariationRealigner.filterSV) and must run before realignment; findDELdisc emits the <DEL>
// breakpoint variations (StructuralVariantsProcessor.findDELdisc) and runs after realignment.
void filterSVStructures(VariationData& vd, int maxReadLength);
// StructuralVariantsProcessor.findDEL: split-read-confirmed <DEL> SVs from the discordant DEL clusters'
// dominant soft clip (findAllSVs runs this FIRST, before findINV/findsv/findDELdisc). Matching the
// soft-clip consensus locates the exact reciprocal breakpoint; the deletion's coverage is then raised
// to the far-breakpoint coverage (via `reload`) so a low-VAF SV on a high-coverage locus is AF-filtered.
// Marks its cluster used so findDELdisc does not re-emit it as a discordant-only estimate.
void findDEL(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
             int maxReadLength, const SVReloadFn& reload);
void findDELdisc(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);

// StructuralVariantsProcessor.findINV: pair-assisted <INV> caller over the discordant same-orientation
// INV clusters. Runs after realignment, BEFORE findsv (so the split-read path skips folded soft clips).
void findINV(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
             int maxReadLength, const SVReloadFn& reload, const std::vector<BamReader*>& bams);
// StructuralVariantsProcessor.findsv: split-read SVs from soft clips. Only the candidate-inversion
// path is emitted here (<INV>); runs after realignment, before findDELdisc (findAllSVs order).
void findsv(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength);
// `reload` re-reads coverage at a large-insertion breakpoint landing outside the region (bi>region.end
// for a 5' clip, bi<region.start for a 3' clip), mirroring the partialPipeline reload so a low-VAF
// duplication/large insertion is AF-filtered instead of pinned at AF=1.0.
void realignlgins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
                  int maxReadLength, const SVReloadFn& reload, const std::vector<BamReader*>& bams);

} // namespace vardict
