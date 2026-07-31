#!/bin/bash
D=~/src/vardictcpp/bench/somatic_dev; cd "$D"
~/src/vardictcpp/build/vardictcpp -G ref.fa -f 0.01 -N "tumor|normal" -b "tumor.bam|normal.bam" \
  -c 1 -S 2 -E 3 -g 4 regions.bed 2>/dev/null | sort > cpp_som.tsv
g=$(wc -l < golden_som.tsv); c=$(wc -l < cpp_som.tsv)
bi=$(comm -12 golden_som.tsv cpp_som.tsv | wc -l)
echo "golden=$g cpp=$c byte-identical-rows=$bi"
echo "--- type counts (col50) golden vs cpp ---"
paste <(cut -f50 golden_som.tsv|sort|uniq -c) <(cut -f50 cpp_som.tsv|sort|uniq -c) 2>/dev/null
