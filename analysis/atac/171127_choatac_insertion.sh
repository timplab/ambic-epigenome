#!/bin/bash

bedpath=/dilithium/Data/Nanopore/Analysis/171025_cho/insertion/bed/insertions.bed
bamdir=/dilithium/Data/NGS/Aligned/171025_choatac/bam
outpath=/dilithium/Data/NGS/Aligned/171025_choatac/insertion/insertion.cov.txt

bampaths=`find $bamdir -name "*nodup.bam"`
header="chr start stop"
for bam in $bampaths;do
  base=$(basename "$bam")
  samp=${base%%.*}
  header="${header} $samp"
done
echo $header | tr " " "\t" > $outpath
samtools bedcov $bedpath $bampaths >> $outpath
