#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/qc

qc=/home/isac/Code/ilee/qc/ontQC.py
qcplotter=/home/isac/Code/ilee/qc/ontQC.R
outdir=/home/isac/Dropbox/Data/ambic/qc
plotpre=$outdir/cho
qcout=$outdir/ontqc.txt
#[ -e $qcout ]&&rm $qcout

for samp in NIHhost NIHIgG SigmaHost SigmaStableGlut SigmaUnstableGlut;do
  for rep in 1 2 3;do
    sum=$(find $root -name "*$samp*rep$rep*")
    echo $samp,$rep
#    out=$(cat $sum | python $qc)
#    printf "%s\t%s\t%s\n" $samp $rep "$out" >> $qcout
  done
done

Rscript $qcplotter -o $plotpre $qcout
