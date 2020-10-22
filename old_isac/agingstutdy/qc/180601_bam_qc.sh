#!/bin/bash
root=/shared/data/analysis
fqdir=$root/fastq
bamdir=$root/pooled/bam
qcroot=$root/qc
logdir=$root/log/qc
[ -e $logdir ]||mkdir -p $logdir
day=0

batchscr="/shared/Code/ilee/slurm/batchcommand.scr"

for samp in Host StableGlut UnstableGlut;do
  for rep in 1 2 3;do
    label=${samp}D${day}rep$rep
    echo $label
    qcdir=$qcroot/$label
    [ -e $qcdir ]||mkdir -p $qcdir
    log=$logdir/$label.qc.log
    sbarg="-c 36 -o $log -e $log $batchscr"
    bam=`readlink -f $bamdir/$label.sorted.bam`
    subcom="NanoPlot -t 36 --bam $bam -o $qcdir -p $label \
      --verbose --raw --store  --drop_outliers \
      -f pdf --N50"
    com="sbatch $sbarg $subcom"
    echo $com
    eval $com
  done
done

