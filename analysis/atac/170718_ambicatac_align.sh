#!/bin/bash

## I'll make this a general wrapper for SE atac-seq

srcdir=/home/isac/Code/timp_genetics/atac
rawdir=/mithril/Data/NGS/Raw/170718_ambicatac/wtimp_choNIHIgGATAC/FASTQ
btidx=/mithril/Data/NGS/Reference/cho/picr/
outdir=/atium/Data/NGS/Aligned/170718_ambicatac
fqdir=${outdir}/fq
subdir=${outdir}/sub
bamdir=${outdir}/bam

if [ ! -d $fqdir ];then
  mkdir $fqdir
  ## for trimming, we have to tell them whether sample is PE or SE,
  ## but just need to provide the input dir
  ${srcdir}/atacTrim.sh -o $fqdir -i ${rawdir}
else
  rm ${outdir}/align.log
  for samp in `readlink -f ${fqdir}/*fq.gz`;do
    echo $samp
    ## for alignment, need to provide second read, but figures out PE by itself
    ${srcdir}/atacAlign.sh -o $bamdir -x $btidx -t 10 $samp &>> ${outdir}/align.log 
  done
f
