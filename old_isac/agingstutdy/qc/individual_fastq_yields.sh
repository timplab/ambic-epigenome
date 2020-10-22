#!/bin/bash
root=/mnt
date=190320
qc=$root/c/qc/${date}_fastq_qc.csv
echo "sample,numreads,numbase,avglen,median,n50" > $qc
fqs=$(find $root/*/reads/* -maxdepth 0 -type f -name "*fastq.gz")
for fq in ; do
  base=$(basename "$fq")
  base=${base%%.*}
#  echo $base
#  gunzip -c $fq | paste - - - - | cut -f2 | wc |\
#    awk -v OFS=',' "{ print \"$base\",\$1,\$3 }" >> $qc
  
done
com="gunzip -c {} | paste - - - - | cut -f2 | wc |\
    awk -v OFS=',' '{ print \"{}\",\$1,\$3 }' >> $qc"
com="python -u /home/ubuntu/Code/ilee/qc/fastq_qc.py \
  {} >> $qc"
parallel $com ::: $fqs
