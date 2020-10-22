#!/bin/bash
root=/shared/data/analysis
pooldir=$root/pooled
bamdir=$root/bam
logroot=$root/log/pool
day=0

batchscr="/shared/Code/ilee/slurm/batchcommand.scr"

for samp in Host StableGlut UnstableGlut;do
  label=${samp}D${day}
  echo $label
  if [ "$1" == "bam" ];then
    outdir=$pooldir/bam
    [ -e $outdir ]||mkdir -p $outdir
    logdir=$logroot/bam
    [ -e $logdir ]||mkdir -p $logdir
    outbam=$outdir/$label.pooled.bam
    bams=`find $outdir -name "*$label*sorted.bam"`
    bamnum=`echo $bams| wc -w`
    if [ "$bamnum" == 1 ];then
      echo "only one bam, just copying"
      com="cp $bams $outbam;cp $bams.bai $outbam.bai"
      echo $com
    else 
      log=$logdir/$label.bampool.log
      subcom="sambamba merge -t 36 -p $outbam $bams"
      sbarg="-c 36 -o $log -e $log $batchscr"
      com="sbatch $sbarg $subcom"
    fi
    echo $com
    eval $com
  fi
done
