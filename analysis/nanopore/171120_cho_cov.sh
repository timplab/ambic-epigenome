#!/bin/bash

bamdir=/dilithium/Data/Nanopore/Analysis/171025_cho/bam/ngmlr
for bam in `find $bamdir -name "*bam"`;do
  base=$(basename "$bam")
  samp=${base%%.*}
  echo $samp
  outname=$samp.cov.txt
  outpath=$bamdir/$outname
  chr=picr_49
  samtools depth -r $chr $bam > $outpath
done
