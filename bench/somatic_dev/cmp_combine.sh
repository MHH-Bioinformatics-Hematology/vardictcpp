#!/bin/bash
# combineAnalysis (merged-BAM refinement) fixture: a tumor-only 15bp deletion (3 explicit-CIGAR reads,
# good var, low coverage) whose supporting reads exist in the normal only as 3' softclips (not promoted
# alone, since realigndel needs a dels5 seed). The merged bam1+bam2 pipeline recovers the normal softclips
# (merged count 7 vs tumor 3), so combineAnalysis reclassifies StrongSomatic -> Germline and back-fills
# the normal block, and FALSE-suppresses the sibling LOH row. golden_comb.tsv is the Java 1.8.3 output.
D=~/src/vardictcpp/bench/somatic_dev; cd "$D"
~/src/vardictcpp/build/vardictcpp -G ref.fa -f 0.01 -N "tumor|normal" -b "tumor2.bam|normal2.bam" \
  -c 1 -S 2 -E 3 -g 4 regions2.bed 2>/dev/null | sort > cpp_comb.tsv
if diff -q golden_comb.tsv cpp_comb.tsv >/dev/null; then
  echo "combineAnalysis fixture: BYTE-IDENTICAL ($(wc -l < golden_comb.tsv) row(s))"
else
  echo "combineAnalysis fixture: DIFFERS"; diff golden_comb.tsv cpp_comb.tsv
fi
