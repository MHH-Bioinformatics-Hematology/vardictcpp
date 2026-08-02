#!/bin/bash
# Self-validating replicate collector for the performance box plots.
# Waits for a quiet node, then times each config; keeps a run ONLY if it is within
# THRESH x its known clean baseline (else it was contended by another job -> discard
# and retry). Accumulates TARGET clean replicates per config. Robust to the shared
# node being intermittently saturated by another benchmark.
set -u
DEV=/home/mh-hannover.local/wolffjoa/src/vardictcpp/bench/wes_parity
V=/home/mh-hannover.local/wolffjoa/src/vardictcpp/build/vardictcpp
JAVA=/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/bin/java
L=/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib
JCP="$L/VarDict-1.8.3.jar:$L/commons-cli-1.2.jar:$L/commons-math3-3.6.1.jar:$L/htsjdk-2.21.1.jar:$L/jregex-1.2_01.jar"
REF=$DEV/ref.fa; T=/tmp/rbc
TARGET=5; THRESH=130   # keep run if wall_s*100 <= THRESH * baseline (i.e. <=1.30x)
MAXTRY=20
RES=$DEV/repbench.tsv; echo -e "sample\tregions\ttool\tthreads\trep\twall_s\tpeakRSS_MB" > $RES

# clean baselines (s) from the clean single-run finalbench: sample:tool:threads -> wall
declare -A BASE=(
 [SRR15006376:cpp:1]=48.5 [SRR15006376:cpp:8]=7.7 [SRR15006376:java:8]=90.3 [SRR15006376:java:1]=185.8
 [SRR15006540:cpp:1]=61.3 [SRR15006540:cpp:8]=9.2 [SRR15006540:java:8]=99.8 [SRR15006540:java:1]=245.2
 [SRR8657348:cpp:1]=133.4 [SRR8657348:cpp:8]=14.3 [SRR8657348:java:8]=142.9 [SRR8657348:java:1]=585.2 )
declare -A REG=( [SRR15006376]=105949 [SRR15006540]=135145 [SRR8657348]=773667 )

met(){ awk -F': ' '/wall clock/{split($2,a,":"); w=(length(a)==3)?a[1]*3600+a[2]*60+a[3]:a[1]*60+a[2]} /Maximum resident/{r=$2/1024} END{printf "%s %s",w,r}' $T; }

wait_quiet(){ # block until no foreign java and load1 < 3
  while :; do
    fj=$(pgrep -x java | wc -l)
    l1=$(awk '{print int($1)}' /proc/loadavg)
    [ "$fj" -eq 0 ] && [ "$l1" -lt 3 ] && { sleep 5; # settle, re-confirm
       fj=$(pgrep -x java|wc -l); [ "$fj" -eq 0 ] && return 0; }
    echo "$(date +%T) waiting for quiet node (foreign_java=$fj load1=$l1)"
    sleep 30
  done
}

runone(){ # sample tool threads cpuset -> echoes "wall rss" or empty if contended
  local s=$1 tool=$2 th=$3 cs=$4
  if [ "$tool" = cpp ]; then
    /usr/bin/time -v taskset -c $cs $V -G $REF -f 0.01 -N $s -b $DEV/$s.bam -c 1 -S 2 -E 3 -g 4 $DEV/$s.bed -th $th >/dev/null 2>$T
  else
    /usr/bin/time -v taskset -c $cs $JAVA -Xmx8g -cp "$JCP" com.astrazeneca.vardict.Main -G $REF -f 0.01 -N $s -b $DEV/$s.bam -c 1 -S 2 -E 3 -g 4 $DEV/$s.bed -th $th >/dev/null 2>$T
  fi
  echo "$(met)"
}

for s in SRR15006376 SRR15006540 SRR8657348; do
  for spec in cpp:1:0 cpp:8:0-7 java:8:0-7 java:1:0; do
    IFS=: read tool th cs <<< "$spec"
    base=${BASE[$s:$tool:$th]}; kept=0; try=0
    while [ $kept -lt $TARGET ] && [ $try -lt $MAXTRY ]; do
      try=$((try+1)); wait_quiet
      read wall rss <<< "$(runone $s $tool $th $cs)"
      # keep if wall*100 <= THRESH*base
      ok=$(awk -v w="$wall" -v b="$base" -v t="$THRESH" 'BEGIN{print (w*100 <= t*b)?1:0}')
      if [ "$ok" = 1 ]; then
        kept=$((kept+1))
        echo -e "$s\t${REG[$s]}\t$([ $tool = cpp ] && echo vardictcpp || echo VarDictJava)\t$th\t$kept\t$wall\t$rss" >> $RES
        echo "$(date +%T) KEEP $s $tool th$th r$kept wall=$wall (base=$base)"
      else
        echo "$(date +%T) DISCARD $s $tool th$th wall=$wall > ${THRESH}%*$base (contended), retry $try"
      fi
    done
  done
done
echo "REPBENCH DONE $(date +%T)"
