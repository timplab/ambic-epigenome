#!/bin/bash
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/summary
qc=../../../qc/ontQC.py
qcplotter=../../../qc/ontQC.R
outdir=/home/isac/Dropbox/Data/ambic/aging_study/qc
qcout=$outdir/181114_ontqc.txt
#[ -e $qcout ]&&rm $qcout

sums=$(find $root -name "*Host1*summary.txt.gz")
for sum in $sums;do
  base=$(basename $sum)
  base=${base%%.*}
  echo $base
#  out=$(gunzip -c $sum | python $qc)
#  printf "%s\t%s\n" $base "$out" >> $qcout
done

qcout=$outdir/ontqc.txt
plotpre=$outdir/181114_chosigma
Rscript $qcplotter -o $plotpre $qcout
