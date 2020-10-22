#!/bin/bash
root=/mnt
rdir=$root/reads
qc=$root/yields_qc.txt
echo "sample readnum basenum" > $qc

for sum in $(find $rdir -maxdepth 1 -name "*summary.txt"); do
  base=$(basename "$sum")
  base=${base%%.*}
  nums=$(awk '{ sum+=$12 }END{ print NR,sum }' $sum)
  echo $base $nums >> $qc
done
