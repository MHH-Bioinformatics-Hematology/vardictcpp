#!/bin/bash
# Benchmark: peak RSS + wall-clock for stock-Java / optimized-Java / vardictcpp on the same
# BAM+reference+region, at -f 0.01 and -f 0. Prints a markdown table.
#
# Usage: run_bench.sh <ref.fa> <in.bam> <region chr:s-e> [chunkSize]
set -u
REF=${1:?ref.fa}; BAM=${2:?in.bam}; REGION=${3:?chr:start-end}; CHUNK=${4:-50000}

L=$HOME/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib
DEPS="$L/commons-cli-1.2.jar:$L/commons-math3-3.6.1.jar:$L/htsjdk-2.21.1.jar:$L/jregex-1.2_01.jar"
STOCK_JAR=${STOCK_JAR:-$L/VarDict-1.8.3.jar}
MOD_JAR=${MOD_JAR:-$L/VarDict-1.8.3.jar}   # override with the optimized jar
VC=${VC:-$(dirname "$0")/../build/vardictcpp}
MAIN=com.astrazeneca.vardict.Main
OUT=$(mktemp -d)

measure() { # label cmd...
  local label="$1"; shift
  /usr/bin/time -v "$@" >"$OUT/o.tsv" 2>"$OUT/t.time"
  local rc=$? rss wall lines
  rss=$(awk '/Maximum resident/{printf "%.2f",$6/1048576}' "$OUT/t.time")
  wall=$(awk -F': ' '/Elapsed/{print $2}' "$OUT/t.time")
  lines=$(wc -l <"$OUT/o.tsv")
  printf "| %-34s | %8s | %8s | %8s |\n" "$label" "${rss}GB" "$wall" "$lines"
}

echo "Region: $REGION   BAM: $(basename "$BAM")"
echo
printf "| %-34s | %8s | %8s | %8s |\n" "tool / setting" "peakRSS" "wall" "lines"
printf "|%s|%s|%s|%s|\n" "------------------------------------" "----------" "----------" "----------"
for f in 0.01 0; do
  echo "**-f $f**"
  measure "stock-Java -Xmx90g"          java -Xmx90g -cp "$STOCK_JAR:$DEPS" $MAIN -G "$REF" -f $f -N s -b "$BAM" -R "$REGION" -th 1
  measure "stock-Java -Xmx8g G1"        java -Xmx8g -XX:+UseG1GC -cp "$STOCK_JAR:$DEPS" $MAIN -G "$REF" -f $f -N s -b "$BAM" -R "$REGION" -th 1
  measure "opt-Java -Xmx8g G1"          java -Xmx8g -XX:+UseG1GC -cp "$MOD_JAR:$DEPS" $MAIN -G "$REF" -f $f -N s -b "$BAM" -R "$REGION" -th 1
  measure "opt-Java -Xmx8g G1 --chunk"  java -Xmx8g -XX:+UseG1GC -cp "$MOD_JAR:$DEPS" $MAIN -G "$REF" -f $f -N s -b "$BAM" -R "$REGION" -th 1 --chunk $CHUNK
  measure "vardictcpp"                  "$VC" -G "$REF" -b "$BAM" -N s -R "$REGION" -f $f
  measure "vardictcpp --chunk"          "$VC" -G "$REF" -b "$BAM" -N s -R "$REGION" -f $f --chunk $CHUNK
done
rm -rf "$OUT"
