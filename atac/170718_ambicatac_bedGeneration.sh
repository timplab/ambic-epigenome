#!/bin/bash

## I'll make this a general wrapper for SE atac-seq

srcdir=/home/isac/Code/timp_genetics/atac
outdir=/atium/Data/NGS/Aligned/170718_ambicatac
nodupout=${outdir}/nodup
bedout=${outdir}/bed

mkdir $bedout
#rm $bedout/*

name="atacNIHIgG"
tag="HNMCFBCXY"

nodups=`readlink -f ${nodupout}/*$tag*.nodup.bam`
${srcdir}/atacbedGeneration.sh -o $bedout -s $name -t $tag $nodups
