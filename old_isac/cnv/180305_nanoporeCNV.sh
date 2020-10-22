#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
bamdir=$root/bam/pooled
beddir=$root/bed
ref=/mithril/Data/NGS/Reference/cho/picr_IgG2/picr_IgG2.fa
refinfo=$ref.chrominfo.txt
outdir=~/Dropbox/Data/Nanopore/180305_cho

if [ "$1" == "genomeprep" ];then
  samtools faidx $ref
  cut -f1,2 $ref.fai > $refinfo
fi

if [ "$1" == "makebed" ];then
  for bam in `find $bamdir -name "*bam" -type f`;do
    base=$(basename "$bam")
    base=${base%.*}
    bedtools bamtobed -i $bam > $beddir/$base.bed
  done
fi


if [ "$1" == "getcov" ];then
  for bed in `find $beddir -name "*host*bed" -type f`;do
    echo $bed
    base=${bed%.bed}
    bedtools genomecov -bga -i $bed -g $refinfo > $base.genomecov.bedGraph
  done
fi

binscript=/home/isac/Code/ilee/sv/binCoverage.py
if [ "$1" == "bincov" ];then
  for bed in `find $beddir -name "*genomecov.bedGraph" -type f`;do
    base=${bed%.bedGraph}
    echo $base
    $binscript -i $bed -w 5000000 -t -1 > $base.binned.5Gb.bedGraph
  done
fi
plotscript=/home/isac/Code/ilee/plot/coverage_plots.R
if [ "$1" == "plot" ];then
  for bincov in `find $beddir -name "*binned.5Gb*" -type f`;do
    echo "plotting"
    base=$(basename "$bincov")
    base=${base%%.bedGraph}
    sampdir=$outdir/$base
    [ -e $sampdir ]||mkdir $sampdir
    out=$sampdir/$base
    echo $base
    Rscript $plotscript genomeCoverage -i $bincov -o $out
  done
fi

