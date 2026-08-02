#!/bin/bash
# Reproduce the GIAB HG002 chr20 exome accuracy benchmark (Table 3, Fig 3).
# Compares VarDict calls (vardictcpp and VarDictJava, which are byte-identical) to the
# NIST v4.2.1 truth with rtg vcfeval, and emits per-region precision/recall/F1.
#
# Requires: samtools >= 1.10 (remote ##idx## support), bedtools, bcftools, rtg-tools,
#           bgzip/tabix, perl (var2vcf_valid.pl from VarDict), python3, and a local hg19.
# Set these to your paths:
set -euo pipefail
HG19=/path/to/hg19.fa                       # chr-prefixed hg19 (chr1..chrY)
VCPP=/path/to/vardictcpp                    # vardictcpp binary
VJAR_CP="VarDict-1.8.3.jar:commons-cli-1.2.jar:commons-math3-3.6.1.jar:htsjdk-2.21.1.jar:jregex-1.2_01.jar"
V2V="perl var2vcf_valid.pl -N HG002 -E -f 0.01"
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
EXO="$GIAB/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome"
EBASE=151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup
TRU="$GIAB/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37"

# 1. chr20 reference named "20" (hg19 primary assembly == GRCh37 sequence for autosomes)
samtools faidx "$HG19" chr20 | sed '1s/^>chr20.*/>20/' > ref20.fa && samtools faidx ref20.fa

# 2. HG002 exome: fetch chr20 reads remotely (index is *.bai, bound via ##idx##)
curl -s "$EXO/$EBASE.bai" -o exome.bai
samtools view -b "$EXO/$EBASE.bam##idx##$PWD/exome.bai" 20 > hg002.exome.20.bam
samtools index hg002.exome.20.bam

# 3. Truth (chr20) + confident BED (chr20)
tabix -h "$TRU/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" 20 > truth20.raw.vcf
curl -s "$TRU/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed" | awk '$1=="20"' > conf20.bed

# 4. Callable target = exome depth>=20 (merged) intersect confident; add a gene column for VarDict
bedtools genomecov -bga -ibam hg002.exome.20.bam | awk '$4>=20' | bedtools merge -d 10 > cov20.bed
bedtools intersect -a cov20.bed -b conf20.bed > target20.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"r"NR}' target20.bed > target20.g.bed

# 5. Call with both implementations; convert to VCF (strandbias.py replaces teststrandbias.R)
for impl in cpp java; do
  if [ "$impl" = cpp ]; then
    "$VCPP" -G ref20.fa -f 0.01 -N HG002 -b hg002.exome.20.bam -c 1 -S 2 -E 3 -g 4 target20.g.bed -th 8
  else
    java -Xmx8g -cp "$VJAR_CP" com.astrazeneca.vardict.Main -G ref20.fa -f 0.01 -N HG002 \
         -b hg002.exome.20.bam -c 1 -S 2 -E 3 -g 4 target20.g.bed -th 8
  fi | sort -k1,1 -k2,2n | python3 strandbias.py | $V2V > $impl.vcf
done
# vardictcpp and VarDictJava VCFs are byte-identical here (verify: diff cpp.vcf java.vcf)

# 6. Normalise, germline operating point (PASS, AF>=0.2), rtg vcfeval over target regions
rtg format -o ref20.sdf ref20.fa
bcftools norm -f ref20.fa -m -any truth20.raw.vcf | bcftools sort | bgzip > truth20.norm.vcf.gz && tabix -fp vcf truth20.norm.vcf.gz
for impl in cpp java; do
  bcftools reheader --fai ref20.fa.fai $impl.vcf | bcftools norm -f ref20.fa -m -any | bcftools sort | bgzip > $impl.norm.vcf.gz
  tabix -fp vcf $impl.norm.vcf.gz
  bcftools view -f PASS $impl.norm.vcf.gz | bcftools view -e 'FMT/AF<0.2' -Oz -o $impl.germ.vcf.gz && tabix -fp vcf $impl.germ.vcf.gz
  rtg vcfeval -b truth20.norm.vcf.gz -c $impl.germ.vcf.gz -t ref20.sdf \
      --evaluation-regions target20.bed -o eval_${impl}_germ --vcf-score-field=QUAL
done

# 7. Per-region distributions for the box plots (writes per_region.tsv)
python3 per_region_metrics.py > per_region.tsv
