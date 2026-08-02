#!/bin/bash
DEV=/home/mh-hannover.local/wolffjoa/src/vardictcpp/bench/wes_parity; V=/home/mh-hannover.local/wolffjoa/src/vardictcpp/build/vardictcpp; JAVA=/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/bin/java; JAR="/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib/VarDict-1.8.3.jar"; DEPS="/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib/commons-cli-1.2.jar:/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib/commons-math3-3.6.1.jar:/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib/htsjdk-2.21.1.jar:/home/mh-hannover.local/wolffjoa/miniconda3/envs/vardict/share/vardict-java-1.8.3-0/lib/jregex-1.2_01.jar"; REF=/home/mh-hannover.local/wolffjoa/src/vardictcpp/bench/wes_parity/ref.fa; T=/tmp/fb
RES=$DEV/finalbench.tsv; echo -e "sample\tregions\ttool\tthreads\twall_s\tpeakRSS_MB" > $RES
met(){ awk -F': ' '/wall clock/{split($2,a,":"); w=(length(a)==3)?a[1]*3600+a[2]*60+a[3]:a[1]*60+a[2]} /Maximum resident/{r=$2/1024} END{printf "%.1f\t%.0f",w,r}' $T; }
runc(){ /usr/bin/time -v taskset -c $3 $V -G $REF -f 0.01 -N $1 -b $DEV/$1.bam -c 1 -S 2 -E 3 -g 4 $DEV/$1.bed -th $2 >/dev/null 2>$T; echo -e "$1\t$4\tvardictcpp\t$2\t$(met)" >>$RES; echo "$(date +%T) cpp $1 th$2 -> $(tail -1 $RES|cut -f5-6)"; }
runj(){ /usr/bin/time -v taskset -c $3 $JAVA -Xmx8g -cp "$JAR:$DEPS" com.astrazeneca.vardict.Main -G $REF -f 0.01 -N $1 -b $DEV/$1.bam -c 1 -S 2 -E 3 -g 4 $DEV/$1.bed -th $2 >/dev/null 2>$T; echo -e "$1\t$4\tVarDictJava\t$2\t$(met)" >>$RES; echo "$(date +%T) java $1 th$2 -> $(tail -1 $RES|cut -f5-6)"; }
for pair in "SRR15006376 105949" "SRR15006540 135145" "SRR8657348 773667"; do set -- $pair; A=$1; R=$2
  runc $A 1 0 $R; runc $A 8 0-7 $R; runj $A 8 0-7 $R; runj $A 1 0 $R
done
echo "FINALBENCH DONE $(date +%T)"
