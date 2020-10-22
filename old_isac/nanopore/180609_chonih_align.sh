#!/bin/bash
env=aws
root=/shared/data/analysis
align=/shared/Code/ilee/oxford/slurm/alignWrapper.sh
batch=/shared/Code/ilee/slurm/batchcommand.scr
ref=/shared/data/Reference/Cricetulus_griseus_chok1gshd.nihIgG.fa
fqdir=$root/fastq
bamdir=$root/bam
logroot=$root/log
pooldir=$root/pooled/bam
[ -e $pooldir ]||mkdir -p $pooldir

if [ "$1" == "align" ];then
  for lab in host IgG;do
    for rep in 1 2 3 ;do
      fq=$(find $fqdir -name "*$lab*$rep.*f*q*gz")
      echo $fq
      base=$(basename "$fq")
      samp=${base%%.*}
      $align -e $env -d $root -b $samp -a ngmlr -r $ref
    done
  done
fi
if [ "$1" == "merge" ];then
  logdir=$logroot/pool
  [ -e $logdir ]||mkdir $logdir
  for samp in IgG host;do
    bams=$(find $bamdir -name "*$samp*.sorted.bam")
    outbam=$pooldir/choNIH$samp.pooled.bam
    subcom="sambamba merge -t 36 -p $outbam $bams"
    log=$logdir/$samp.pool.log
    sbatch -c 36 -e $log -o $log $batch $subcom
  done
fi
