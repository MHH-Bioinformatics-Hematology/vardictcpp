# Known VarDictJava parity divergences (vardictcpp)

Verified against Java `-th 1` golden output. Each entry records a row where cpp does not
(yet) match VarDictJava, the read-level root cause, and why it is still open.

The reference dataset is the CCLE WES sample **SRR8657348** (hg19, full BED), where cpp is
byte-identical to the Java `-th 1` golden on **54448 / 54448** rows. The two previously-open
structural-variant divergences (a spurious split-read `<DUP>` at chr18:18518111 and a missing
discordant `<DUP>` at chr20:58889503) are now **closed**; see the CLOSED section below.

There are currently **no known open divergences** on the verified WES parity samples
(SRR8657348 54448/54448; SRR15006386, SRR15006375, SRR15006376, SRR15006540 all byte-identical
to Java `-th 1`).

---

## CLOSED — chr18:18518111 spurious split-read `<DUP>` (+ co-located SNV `Seg` pollution)

Fixed in two coupled commits:

1. **Reference-map extent** (`reference: drop the trailing SEED_1 of a non-end window from
   has()`). `ReferenceResource.getReference` fills `referenceSequences` only up to
   `sequenceEnd - SEED_1` when the loaded window does not reach the contig end
   (`ReferenceResource.java:106`), so `ref.get(p) == null` for the trailing `SEED_1 (=17)`
   bases. `Reference::has()` now returns false past `loadedEnd() - SEED_1` (tracked as
   `primaryEffEnd_`) unless the window reaches the contig end; every `at()` consumer guards with
   `has()`, so `ismatchref` rejects the soft-clip seed match at 18519646 exactly as Java does and
   no `getSV(...,DUP)` is fabricated. This also clears the `Seg 10-0-0` pollution on the
   co-located `A>G` SNV.

2. **Inter-chromosomal fusion clusters + SOFTP2SV guard** (`realigner: port inter-chr fusion
   clusters + SOFTP2SV findsv guard`). Once the `<DUP>` is gone, the 18518112 soft clip is no
   longer `used`, so `findsv` fabricated a split-read `<INV>` (chr18:18517619-18518122). Java
   suppresses it via `SOFTP2SV` (`StructuralVariantsProcessor.java:705`): the soft clip is
   claimed by a discordant SV cluster whose `checkCluster` came out degenerate (`used=true`).
   That cluster is an inter-chromosomal fusion cluster (`svrfus["chr12"]`, 4 chr12 mates) built
   in CigarParser's inter-chr branch (`CigarParser.java:2205-2269`), which vardictcpp skipped
   (early return on `mtid != tid`). It is now ported (svffus/svrfus keyed by mate contig, with
   the nearby DEL/DUP/INV disc bumps), every SV cluster's dominant soft clip registers in
   `SOFTP2SV` during `filterSV`, and both `findsv` soft-clip loops honour the
   `SOFTP2SV[p][0].used` guard. cpp now builds the `cnt=4 softp=18518112 used=true` cluster
   exactly like Java and drops the spurious `<INV>`.

## CLOSED — chr20:58889503 discordant+split `<DUP>`

Fixed by `realigner: port findDUPdisc (discordant-pair-supported <DUP>)`.
`StructuralVariantsProcessor.findDUPdisc` (`StructuralVariantsProcessor.java:1330-1608`) folds
discordant read-pair support into a `<DUP>` at a duplication cluster's breakpoint: it refines the
breakpoint from the dominant soft clip's split-read match, builds the `+<dup>` insertion, and
sums the pair count plus both soft-clip split counts into the SV marker. vardictcpp collected
svfdup/svrdup clusters for disc bookkeeping only and never consumed them. Both the forward
(svfdup) and reverse (svrdup) paths are now ported (mirroring the realignlgins split-read DUP
emission and findDELdisc's cluster gating) and run last in `findAllSVs`. This recovers
`chr20 58889503 58889880 C <DUP>` (Seg `34-10-2`, 10 discordant pairs).
