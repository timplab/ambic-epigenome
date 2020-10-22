#!/bin/bash
day=0
root=/shared/data/analysis
pooldir=$root/pooled
bamdir=$pooldir/bam
beddir=$pooldir/bed
[ -e $beddir ]||mkdir $beddir
batch=/shared/Code/ilee/slurm/batchcommand.scr
logdir=$root/log/bed
[ -e $logdir ]||mkdir $logdir

for bam in `find $bamdir -name "*NIH*bam" -type f`;do
  fn=$(basename "$bam")
  prefix=${fn%.bam*}
  bed=$beddir/$prefix.bed.gz
  echo "making bed file for $prefix"
  com="bedtools bamtobed -i $bam |\
    sort -T $beddir -k1,1 -k2,2n |\
    bgzip > $bed;\
    tabix -p bed $bed"
  log=$logdir/$prefix.bamtobed.log
  sbatch -c 36 -e $log -o $log $batch $com
  

done
