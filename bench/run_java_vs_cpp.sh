#!/bin/bash
# Head-to-head: UNFIXED stock VarDictJava 1.8.3 vs vardictcpp on real AML-panel data.
# Measures peak RSS + wall-clock (/usr/bin/time -v) at -f 0.01 and -f 0, single- and multi-thread,
# then checks output equivalence (shared variant set + byte-identical rows) at the -f 0.01 calling threshold.
set -u
REF=$HOME/transvar_test/hg19.fa
BAM=/mnt/processed/wolffjoa/mrd_redo/bwa_mem/MRD2020-032_H20PB2025.bwa.bam
BED=$(dirname "$0")/panel4.bed
VC=$(dirname "$0")/../build/vardictcpp
L=$HOME/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib
DEPS="$L/commons-cli-1.2.jar:$L/commons-math3-3.6.1.jar:$L/htsjdk-2.21.1.jar:$L/jregex-1.2_01.jar"
JAR=$L/VarDict-1.8.3.jar    # the UNFIXED stock 1.8.3 jar (no memory/throughput patches)
MAIN=com.astrazeneca.vardict.Main
T=$(mktemp -d)

# UNFIXED Java = default JVM (no -Xmx tuning, no G1, no --chunk): the out-of-the-box experience.
jvar() { java -cp "$JAR:$DEPS" $MAIN -G "$REF" -N s -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$BED" "$@"; }
cvar() { "$VC" -G "$REF" -N s -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$BED" "$@"; }

measure() { # outfile label cmd...
  local out="$1" label="$2"; shift 2
  /usr/bin/time -v "$@" >"$out" 2>"$T/t"
  local rss wall; rss=$(awk '/Maximum resident/{printf "%.3f",$6/1048576}' "$T/t")
  wall=$(awk -F': ' '/wall clock/{print $2}' "$T/t")
  local lines; lines=$(grep -c . "$out")
  printf "| %-24s | %9s | %10s | %7s |\n" "$label" "${rss} GB" "$wall" "$lines"
}

echo "# vardictcpp vs UNFIXED VarDictJava 1.8.3 - real AML panel (MRD2020-032, hg19, 698 regions)"
echo "Host: $(nproc) cores; BAM $(du -h "$BAM"|cut -f1); ref hg19; stock jar $(basename "$JAR")"
echo
for f in 0.01 0; do
  for th in 1 4; do
    echo "## -f $f, -th $th"
    printf "| %-24s | %9s | %10s | %7s |\n" "tool" "peak RSS" "wall" "lines"
    printf "|%s|%s|%s|%s|\n" "--------------------------" "-----------" "------------" "---------"
    measure "$T/j_${f}_${th}.tsv" "stock Java (default JVM)" java -cp "$JAR:$DEPS" $MAIN -G "$REF" -N s -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$BED" -f $f -th $th
    measure "$T/c_${f}_${th}.tsv" "vardictcpp"               "$VC" -G "$REF" -N s -b "$BAM" -c 1 -S 2 -E 3 -g 4 "$BED" -f $f -th $th
    echo
  done
done

echo "## Output equivalence at -f 0.01 (-th 4, the calling threshold)"
J=$T/j_0.01_4.tsv; C=$T/c_0.01_4.tsv
key() { awk -F'\t' '{print $3"_"$4"_"$6"_"$7}' "$1" | sort; }
comm -23 <(key "$J") <(key "$C") > "$T/fn"   # in Java, not C++
comm -13 <(key "$J") <(key "$C") > "$T/fp"   # in C++, not Java
jn=$(grep -c . "$J"); cn=$(grep -c . "$C")
bi=$(comm -12 <(sort "$J") <(sort "$C") | wc -l)
shared=$(comm -12 <(key "$J") <(key "$C") | wc -l)
echo "Java variants: $jn   C++ variants: $cn"
echo "Shared variant keys (chr/pos/ref/alt): $shared"
echo "Full-row byte-identical: $bi / $jn"
echo "Set-level FN (Java-only): $(grep -c . "$T/fn")   FP (C++-only): $(grep -c . "$T/fp")"
[ -s "$T/fn" ] && { echo "  FN keys:"; sed 's/^/    /' "$T/fn"; }
[ -s "$T/fp" ] && { echo "  FP keys:"; sed 's/^/    /' "$T/fp"; }
cp "$J" "$T/../j_panel.tsv" 2>/dev/null
echo
echo "(raw outputs in $T)"
