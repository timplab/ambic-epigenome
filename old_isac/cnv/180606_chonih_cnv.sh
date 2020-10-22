#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
beddir=$root/bed
svdir=$root/sv
bin=/home/isac/Code/ilee/sv/binBedcounts.py
plotdir=~/Dropbox/Data/ambic/plots
cnv=/home/isac/Code/ilee/sv/cnv.R
win=10000
faidx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.nihIgG.fa.fai

beds=`find $beddir -maxdepth 1 -name "*bed.gz" -type f`
for bed in $beds;do
  base=$(basename "$bed")
  echo $base
  pre=${base%%.*}
  counts=$svdir/$pre.counts.txt
  if [ ! -e $counts ];then
    python $bin -i $bed -w $win > $counts
  fi
done
counts=`find $svdir -name "*counts.txt"`
echo $counts
Rscript $cnv -w $win -o $plotdir/choNIH_cnv.pdf -f $faidx $counts

