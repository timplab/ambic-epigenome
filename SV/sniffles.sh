#!/bin/bash
#root=/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/oldgenome
#bam=${root}/170608_choIgGNIH_hits.sorted.bam
bam=`readlink -f $1`
root=$(dirname "$bam")
base=${1%.*}

echo bam:${bam}

~/Code/Sniffles/bin/*/sniffles -m $bam \
			       -s 2 -q 5 -n 10 \
			       -v ${root}/${base}_sniffles.vcf &> ${base}_sniffles.log
