#!/bin/bash

wdir=/shared/insertion
bamdir=$wdir/insertbam
beddir=$wdir/bed
pname=4_0cdhfr_vrc01wtg1m3_dgv

for bam in `find $bamdir -name "*bam" -type f`;do
  dir=$(dirname "$bam")
  base=$(basename "$bam")
  samp=${base%%.*}
  bedout=$beddir/$samp.insertion.bed
  #bedtools bamtobed -i $bam > $bedout
  echo $samp
  awk -v p="$pname" \
    '{if($1!=p && $3-$2>50) print $1,$2,$3 }' $bedout |\
    wc -l
done
