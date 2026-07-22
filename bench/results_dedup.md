# vardictcpp vs VarDictJava 1.8.3 — UMI-deduplicated MRD panel

Real-data validation on **UMI-deduplicated** AML MRD amplicon BAMs (the production dedup output,
`fgbio`/`umi_tools`-collapsed consensus reads), which is the input the caller actually runs on in
the MRD pipeline. Both tools run the **same** command over the same 698-region hg19 panel.

- Data: `/mnt/processed/wolffjoa/mrd_redo/umi_dedup/*.dedup.bam` (hg19, coordinate-sorted), 12 largest samples
- Reference: `hg19.fa`; Panel: 698 covered regions (`panel4.bed`); `-f 0.01`, `-th 8`
- Java: the `optimize-memory` branch jar (byte-identical to stock 1.8.3); C++: `build/vardictcpp`
- Metrics: variant-set FP/FN vs Java (key = chr,pos,ref,alt), full-row byte-identity, Depth/AltDepth
  agreement on shared variants, peak RSS + wall-clock (`/usr/bin/time -v`).

| Sample | Java | C++ | FP | FN | byte-ident | Depth/AltD | Java RSS | C++ RSS | Java wall | C++ wall |
|---|---|---|---|---|---|---|---|---|---|---|
| MRD2019-035_H19KM1798 | 1 | 1 | 0 | 0 | 1/1 | 1/1 | 0.89 GB | 0.024 GB | 2.1 s | 0.5 s |
| MRD2019-054_H19KM2702 | 1 | 1 | 0 | 0 | 1/1 | 1/1 | 0.50 GB | 0.025 GB | 2.3 s | 0.5 s |
| MRD2019-079_H19KM3532 | 0 | 0 | 0 | 0 | — | — | 0.62 GB | 0.025 GB | 2.2 s | 0.4 s |
| MRD2020-060_H18PB1757 | 0 | 0 | 0 | 0 | — | — | 0.90 GB | 0.022 GB | 1.8 s | 0.4 s |
| MRD2019-040_H19KM2056 | 2 | 3 | 1 | 0 | 1/2 | 1/2 | 0.90 GB | 0.024 GB | 2.0 s | 0.5 s |
| MRD2019-037_H19KM1814 | 0 | 0 | 0 | 0 | — | — | 1.23 GB | 0.025 GB | 6.8 s | 0.9 s |
| MRD2020-034_H17KM1413 | 0 | 0 | 0 | 0 | — | — | 0.89 GB | 0.026 GB | 2.2 s | 0.5 s |
| MRD2023-005_H23KM1055 | 0 | 0 | 0 | 0 | — | — | 0.89 GB | 0.021 GB | 2.0 s | 0.4 s |
| MRD2023-014_H23PB1357 | 0 | 0 | 0 | 0 | — | — | 0.45 GB | 0.021 GB | 2.0 s | 0.4 s |
| MRD2025-006_H25KM1467 | 0 | 0 | 0 | 0 | — | — | 0.92 GB | 0.022 GB | 2.5 s | 0.4 s |
| MRD2023-006_H23KM1033 | 0 | 0 | 0 | 0 | — | — | 0.54 GB | 0.025 GB | 2.4 s | 0.4 s |
| MRD2019-071_H19KM3257 | 3 | 4 | 1 | 0 | 2/3 | 2/3 | **2.57 GB** | **0.028 GB** | **230 s** | **14.4 s** |

## Accuracy

- **10 of 12 samples: exact parity** — identical variant set (0 FP / 0 FN) and byte-identical rows.
- **0 false-negatives across all 12 samples.** Every variant VarDict calls, vardictcpp also calls.
- **2 false-positives total**, both on the two samples that carry them, and both are the **single
  documented `CigarModifier` gap** (reads VarDict rewrites-before-counting that this port counts):
  - `chr2:33141508 A>G` (MRD2019-040): 2 poly-G MiSeq artifact reads (`110M81S`, `179M29S`) whose
    G>A tail VarDict soft-clips away via `captureMisSoftly3Mismatches`; this port counts the SNV.
  - `chr15:90631849 GA>G` (MRD2019-071): a spurious 1 bp deletion in `15S9M1D126M` reads — a 9 bp
    mapped island precedes the indel behind a 15 bp soft-clip; VarDict reshapes the read so the
    deletion never counts.

### Deletion representation fixed in this pass

On MRD2019-040 the shared deletion `chr15:90631879 TG>T` previously diverged (C++ emitted a
non-normalized `GG>G`, type `Complex`). Ported VarDict's exact deletion handling
(`ToVarsBuilder` + `proceedVrefIsDeletion`): anchor one base 5′ of the stored position
(`ref[p-1] + deleted bases` / `ref[p-1]`, `startPosition--`), deletion `findMSI`, the `genotype1/-N`
genotype form, `varType`-based `Deletion` classification, and `NM = edit_distance − (I+D length)`.
The row is now **byte-identical to VarDict on 30 of 32 columns**; the only remaining difference is
`Depth` (52 vs 54) — the documented distributed-coverage accounting at indel-dense positions — which
cascades to `AF`. This fix also eliminated **every** prior deletion false-negative (the mismatched
normalization had made the same event look like a different key).

## Performance (the original motivation: `-f 0` → ~200 GB Java)

Even on these small deduplicated panels the memory/latency gap is large, and it widens sharply on the
one sample with a repeat-dense pile-up (MRD2019-071):

- **Peak RSS: 39× lower on average, 91× lower at peak** (2.57 GB → 0.028 GB on MRD2019-071).
- **Wall-clock: 13× faster on average, 16× faster on the slowest sample** (230 s → 14.4 s).
- C++ peak RSS is essentially flat (~24 MB) regardless of pile-up depth, because per-region
  structures are released as soon as the region is emitted — the property that removes the Java
  `-f 0` blow-up entirely.

## Regression guard

The deletion/genotype/NM changes were verified not to regress prior parity:
- Self-contained fixture (`test/run_tests.sh`): **100 % byte-identical** (0 FP/FN).
- 1 Mb / 300× synthetic (`sim:1-1000000`, `-f 0.01`): 1904 vs 1904 variants, **98.0 % byte-identical**,
  Depth/AltDepth exact on 1875/1901; the residual FP/FN are the pre-existing MNV-growing edge cases
  (all multi-base MNVs, e.g. `TGTG>CGTC`), unrelated to the deletion path.
