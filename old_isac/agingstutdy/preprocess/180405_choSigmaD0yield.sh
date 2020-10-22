#!/bin/bash

root=/shared/data/analysis/choSigma
fqroot=$root/fastq

for dir in `find $fqroot/* -maxdepth 0 -type d`;do
  samp=$(basename "$dir")
  yield=`awk 'OFS=","{ sum+=$13 }END{ print sum,NR }' $dir/*summary*`
  echo $samp,$yield

done
