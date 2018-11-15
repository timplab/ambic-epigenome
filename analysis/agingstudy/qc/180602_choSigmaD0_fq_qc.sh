#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
fqdir=$root/fastq
metadata=$root/choSigmaAgingStudy_data.csv
qc=/home/isac/Code/ilee/qc/ontQC.py
qcplotter=/home/isac/Code/ilee/qc/ontQC.R
outdir=/home/isac/Dropbox/Data/ambic/aging_study/qc
plotpre=$outdir/chosigmaD$day
qcout=$outdir/ontqc.txt
#[ -e $qcout ]&&rm $qcout

awk 'NR>1' $metadata | while IFS=$',' read -r -a line
do
  n=$(($n+1)) 
#  if [ $n -gt 1 ];then
#    exit
#  fi
  lab=${line[1]}
  run=${lab#*choSigma}
  samp=${run%D0*}
  rep=${run#*rep}
  rep=${rep%_*}
  sum=`readlink -f $fqdir/$lab*summary.txt`
#  out=$(python $qc $sum)
#  printf "%s\t%s\t%s\n" $samp $rep "$out" >> $qcout
done

Rscript $qcplotter -o $plotpre $qcout
