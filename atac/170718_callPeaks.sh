#!/bin/bash

## I'll make this a general wrapper for SE atac-seq

srcdir=/home/isac/Code/timp_genetics/atac
outdir=/atium/Data/NGS/Aligned/170718_ambicatac
bedout=${outdir}/bed
peakout=${outdir}/peak

mkdir $peakout
rm -r $peakout/*

name="atacNIHIgG"
tag="HNMCFBCXY"

beds=`readlink -f ${bedout}/*$tag*bed.gz`

for bed in $beds; do
  s=${bed##*/}
  s=${s%.bed.gz}
  ${srcdir}/atacPeak.sh -o $peakout -s $s -g mm $bed

done


