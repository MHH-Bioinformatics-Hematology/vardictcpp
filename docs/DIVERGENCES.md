# Known VarDictJava parity divergences (vardictcpp)

Verified against Java `-th 1` golden output. Each entry records a row where cpp does not
(yet) match VarDictJava, the read-level root cause, and why it is still open.

The reference dataset is the CCLE WES sample **SRR8657348** (hg19, full BED), where cpp is
byte-identical to the Java `-th 1` golden on **54446 / 54448** rows. The two remaining
divergences below are both deep structural-variant (SV) creation/ordering issues. Fixing
either faithfully requires a multi-part SV rework whose regression risk (against the many
SV rows that are already byte-identical) is not justified by the two rows it would recover,
so they are documented rather than patched.

Both are **cpp bugs** (Java is correct here), not Java bugs.

---

## OPEN 1 — chr18:18518111 spurious split-read `<DUP>` (+ co-located SNV `Seg` pollution)

**Rows (cpp only, Java has neither the DUP nor the non-zero Seg):**

```
cpp:    chr18 18518111 18518111 A G       Seg 10-0-0     (SNV, should be Seg 0)
cpp:    chr18 18518111 18519645 A <DUP>   Seg 10-0-0     (spurious duplication, A/+1535)
golden: chr18 18518111 18518111 A G       Seg 0
```

**Root cause (two coupled cpp bugs; the second is masked by the first):**

1. **Reference-map extent.** VarDictJava's `ReferenceResource.getReference` stores
   `referenceSequences` only up to `sequenceEnd - SEED_1` when the loaded window does **not**
   reach the chromosome end (`siteEnd = exon.length() - SEED_1`, `ReferenceResource.java:106`).
   The trailing `SEED_1 (=17)` bases of the padded window are therefore absent from the
   reference map (`ref.get(p) == null` there). vardictcpp's `Reference::has()/at()` instead
   expose the **whole** loaded window (`reference.hpp`), so cpp answers `has()` for ~17 bases
   past where Java's map ends.

   For the cnt-10 5' soft clip at 18518112 (consensus reversed
   `CAAAAAGAGTGTTTCCAAACTGCTGCATC`), `realignlgins` → `findMatch` seeds at 18519626 and calls
   `ismatchref(seq, 18519646, dir=-1, MM=1)`. The true hg19 reference at 18519617–18519646 is
   `acaaaaagagtgtttccaaactgctgtatc`, giving exactly **1 mismatch** (ref `t` vs clip `C` at
   18519643) → `mm=1 ≤ 1` → cpp accepts and builds `getSV(...,DUP)` at `bi=p-1=18518111`.
   Java's `ismatchref` hits `ref.get(18519646) == null` (its map ends near 18519632) → returns
   `false` → `findMatch` returns 0 → `realignlgins` `continue`s (`VariationRealigner.java:1645`),
   **no DUP**. The `Seg 10-0-0` on the co-located `A>G` SNV is a side effect: cpp's DUP marker
   at 18518111 makes the SNV inherit the position-level `splits-pairs-clusters` string.

   *A faithful reference-extent truncation (cap `has()/at()` and the seed map at
   `loadedEnd - SEED_1` unless the window reaches the chromosome end) was implemented and
   verified to suppress the DUP AND fix the SNV Seg, while keeping all gates + the four
   byte-identical samples unchanged. It was reverted only because it unmasks bug (2).*

2. **Missing discordant DUP cluster (unmasked once the DUP is gone).** With the DUP no longer
   created, the 18518112 soft clip is no longer marked `used`, so `findsv` reaches it and its
   reverse-strand `findMatchRev` fabricates a split-read `<INV>` at chr18:18517619-18518122
   (`Seg 10-0-0`). Java suppresses this via the `SOFTP2SV` guard
   (`StructuralVariantsProcessor.java:705`: skip a soft clip whose SV cluster is `used`).
   The relevant cluster is Java's discordant `svrdup`/DEL cluster with `Softp: 18518112`,
   `Used: true` (from the discordant read pairs at 18518112, e.g. `8503378 28S40M` MQ 14,
   isize 1844). vardictcpp builds **no** SV cluster at that softp: it drops the MQ-14 mate
   (`prepareSVStructures` `MQ < 15` gate) and does not populate the reverse DUP cluster's
   `soft` map for this pair, so its `SOFTP2SV` equivalent has no entry to gate on.

   *A `SOFTP2SV` map + `findsv` guard were also implemented (faithful to Java), but they only
   help once cpp actually builds the underlying discordant cluster — a parse-level SV-structure
   port that is the larger, riskier piece.*

**Why open:** a correct fix needs BOTH the reference-extent truncation AND porting the missing
discordant DUP-cluster building (+ `SOFTP2SV` gate). The first alone trades the DUP FP for an
INV FP. The parse-level SV-cluster work touches `prepareSVStructures`/`filterSV`, which feed
many already-byte-identical SV rows, so it is deferred as a dedicated SV-structure task.

---

## OPEN 2 — chr20:58889503 discordant+split `<DUP>` missing (1 FN)

**Row (golden only):**

```
golden: chr20 58889503 58889880 C <DUP>  Seg 34-10-2   (C/+378, discordant DUP)
cpp:    (absent)
```

**Root cause:** the `Seg 34-10-2` shows this duplication is carried by discordant read-pair
support: 34 splits, **10 pairs**, 2 clusters. Java emits it through the discordant-DUP fold-in
(`StructuralVariantsProcessor.findDUPdisc` + the pair-supported DUP branch of the split-read
path). vardictcpp deliberately does **not** port `findDUPdisc` (see `cigar_parser.hpp`: "DUP
clusters are collected only for the exact cross-orientation disc bookkeeping; findDUPdisc is
not ported"). Without it there is no path to fold discordant-pair support into a DUP, so the
call is never made.

**Why open:** porting `findDUPdisc` (and the pair-support fold-in that raises the DUP's
depth/AF above threshold) is a self-contained but non-trivial SV feature. It affects the
discordant-DUP path shared with other loci, so it is deferred as a dedicated task rather than
risk the byte-identical SV rows for this single FN.
