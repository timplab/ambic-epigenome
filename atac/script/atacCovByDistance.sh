#!/bin/bash

module=$1
shift 1
while :
do
  case "$1" in
    -i | --indir)
      beddir=$2 # input bed directory
      shift 2
      ;;
    -a | --annotation)
      ann=$2 #input annotation bed
      shift 2
      ;;
    -t | --tag)
      tag=$2
      shift 2
      ;;
    *) break
      ;;
  esac
done

outdir=$beddir/$tag
[ -e $outdir ]||mkdir $outdir

if [ "$module" == "getdist" ];then
  for bed in `find $beddir -name "*bed.gz"`;do
    base=$(basename "$bed")
    base=${base%%.*}
    echo $base
    gunzip -c $bed |\
      bedtools closest -D b -b $ann -a stdin |\
      awk '{ OFS="\t" }{ if(($NF>-2000)&&($NF<2000)) print }'> $outdir/$base.$tag.bed
    # parse
    awk '{ OFS="\t" }{ print $NF,$10 }' $outdir/$base.$tag.bed > $outdir/$base.${tag}dist.txt
  done
fi

if [ "$module" == "plotdist" ];then
  for bed in `find $beddir -name "*bed.gz"`;do
    base=$(basename "$bed")
    base=${base%%.*}
    echo $base
    Rscript ~/Code/ambic-epigenome/plot/atacSeq_plots.R covByDistance -i $outdir/$base.${tag}dist.txt -o $outdir/$base.${tag}dist.pdf
  done
fi
