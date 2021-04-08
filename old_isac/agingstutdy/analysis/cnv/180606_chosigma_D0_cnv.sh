#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
pooldir=$root/pooled
beddir=$pooldir/bed
svdir=$root/sv
bin=/home/isac/Code/ilee/sv/binBedcounts.py
plotdir=~/Dropbox/Data/ambic/aging_study/plots
cnv=/home/isac/Code/ilee/sv/cnv.R
win=10000
faidx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.fa.fai

beds=`find $beddir -name "*bed.gz" -type f`
for bed in $beds;do
  base=$(basename "$bed")
  echo $base
  pre=${base%%.*}
  counts=$svdir/$pre.counts.$win.txt
  if [ ! -e $counts ];then
    python $bin -i $bed -w $win > $counts
  fi
done
counts=`find $svdir -name "*counts.$win.txt"`
echo $count
Rscript $cnv -w $win -o $plotdir/day${day}_cnv.pdf -f $faidx $counts
