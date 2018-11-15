#!/bin/bash
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/summary
qc=../../../qc/ontQC.py
qcplotter=/home/isac/Code/ilee/qc/ontQC.R
outdir=/home/isac/Dropbox/Data/ambic/aging_study/qc
plotpre=$outdir/181114_chosigma
qcout=$outdir/181114_ontqc.txt
[ -e $qcout ]&&rm $qcout

sums=$(find $root -name "*summary.txt.gz")
for sum in $sums;do
  base=$(basename $sum)
  base=${base%%.*}
  out=$(gunzip -c $sum | python $qc)
  echo $base
  printf "%s\t%s\n" $base "$out" >> $qcout
done

exit

Rscript $qcplotter -o $plotpre $qcout
