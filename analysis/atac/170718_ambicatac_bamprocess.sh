#!/bin/bash

## I'll make this a general wrapper for SE atac-seq

srcdir=/home/isac/Code/timp_genetics/atac
outdir=/atium/Data/NGS/Aligned/170718_ambicatac
nodupout=${outdir}/nodup
bamdir=${outdir}/bam

mkdir $nodupout
#rm $nodupout/*

for bam in `ls ${bamdir}/*.sorted.bam`;do
  samp=${bam##*/}
  samp=${samp%%.*}
  ${srcdir}/atacbamProcess.sh -o $nodupout $bam 
  
done
