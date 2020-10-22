#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
pooldir=$root/pooled
bamdir=$pooldir/bam
beddir=$pooldir/bed
[ -e $beddir ]||mkdir $beddir

for bam in `find $bamdir -name "*bam" -type f`;do
  fn=$(basename "$bam")
  prefix=${fn%.bam*}
  bed=$beddir/$prefix.bed.gz
  if [ ! -e $bed ];then
    echo "making bed file for $prefix"
    bedtools bamtobed -i $bam |\
      sort -T $beddir -k1,1 -k2,2n |\
      bgzip > $bed
    tabix -p bed $bed
  fi
done
