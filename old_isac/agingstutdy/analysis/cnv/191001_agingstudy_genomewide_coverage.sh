#!/bin/bash
root=/kyber/Data/Nanopore/projects/ambic/sigma
bdir=$root/bed/pooled_rep
bindir=$root/coverage/pooled_rep

w=10000
bins=$bindir/picr_sigma_bins.10kb.bed
if [ "$1" == "bingenome" ];then
  binner=/home/isac/Code/ilee/util/bin_genome.py
  ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.sigmaIgG.fa
  $binner -g $ref -b $w > $bins
fi

bases=$(find $bdir -name "*bed.gz" | tr "/" " " | awk '{print $NF}' | cut -d"." -f1)

if [ "$1" == "cov" ];then
  [ -e $bindir ]||mkdir $bindir
  base="{}"
  bed=$bdir/$base.pooled.bed.gz
  cov=$bindir/$base.binned_coverage.$w.bed
  com="bedtools coverage -a $bins -b $bed > $cov"
  echo $com
  parallel $com ::: $bases
fi

if [ "$1" == "merge" ];then
  scr=~/Code/ambic-epigenome/scripts/make_coverage_matrix_from_bedtools_coverage.sh
  out=$bindir/CHOZN_coverage_matrix.$w.bed
  for base in $bases; do
    cov=$bindir/$base.binned_coverage.$w.bed
    inputs="$inputs $cov"
  done
  $scr $inputs > $out
fi
  
