#!/bin/bash
root=/atium/Data/Nanopore/Analysis
rawdir=/atium/Data/Nanopore/oxford
srcdir=/home/isac/Code/ilee/oxford

for runpath in `ls -d ${rawdir}/*choIgGNIH*`; do
  base=$(basename "$runpath")
  echo $base
  fastqpath=`ls ${runpath}/${base}*fq.gz`
  summarypath=`ls ${runpath}/${base}*.csv.gz`
  statpath=${root}/${base}/${base}_stats.txt

  ${srcdir}/oxford_stats.R $summarypath $statpath
  
done
