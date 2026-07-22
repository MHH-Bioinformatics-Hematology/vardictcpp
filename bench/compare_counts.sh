#!/bin/bash
# Validate vardictcpp counting against VarDictJava pileup on matched SNV alleles.
# Usage: compare_counts.sh <ref.fa> <in.bam> <region chr:s-e>
set -u
REF=${1:?ref.fa}; BAM=${2:?in.bam}; REGION=${3:?chr:s-e}
L=$HOME/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib
DEPS="$L/commons-cli-1.2.jar:$L/commons-math3-3.6.1.jar:$L/htsjdk-2.21.1.jar:$L/jregex-1.2_01.jar"
JAR=${STOCK_JAR:-$L/VarDict-1.8.3.jar}
VC=${VC:-$(dirname "$0")/../build/vardictcpp}
T=$(mktemp -d)

java -Xmx8g -cp "$JAR:$DEPS" com.astrazeneca.vardict.Main -G "$REF" -f 0.01 -N s -b "$BAM" -R "$REGION" -p -th 1 > "$T/j.tsv" 2>/dev/null
"$VC" -G "$REF" -b "$BAM" -N s -R "$REGION" -f 0.01 -p > "$T/c.tsv" 2>/dev/null

awk -F'\t' 'length($6)==1&&length($7)==1&&$6!=$7{print $4"_"$6"_"$7"\t"$8"\t"$9}' "$T/j.tsv" | sort > "$T/j.snv"
awk -F'\t' 'length($6)==1&&length($7)==1&&$6!=$7{print $4"_"$6"_"$7"\t"$8"\t"$9}' "$T/c.tsv" | sort > "$T/c.snv"
join "$T/j.snv" "$T/c.snv" | awk '{n++; if($2==$4)d++; if($3==$5)a++; dd+=($2>$4?$2-$4:$4-$2); ad+=($3>$5?$3-$5:$5-$3)}
  END{printf "matched SNV alleles: %d\nDepth exact:    %d/%d (%.1f%%)  meanAbsDiff %.3f\nAltDepth exact: %d/%d (%.1f%%)  meanAbsDiff %.3f\n", n, d,n,100*d/n,dd/n, a,n,100*a/n,ad/n}'
rm -rf "$T"
