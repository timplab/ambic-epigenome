#!/bin/bash
batch=/shared/Code/ilee/slurm/batchcommand.scr
root=/shared/data/analysis
logdir=$root/log/sniffles
[ -e $logdir ]||mkdir $logdir
bamdir=$root/pooled/bam
bams=$(find $bamdir -name "*pooled.bam" -type f)
outdir=$root/sniffles
[ -e $outdir ]||mkdir $outdir

for bam in $bams;do
  base=$(basename "$bam")
  label=${base%%.*}
  out=$outdir/$label.sniffles.vcf
  if [ -e $out ];then
    continue
  fi
  echo $label
  log=$logdir/$label.sniffles.log
  com="sniffles -m $bam -s 2 -n 30 \
    -v $out"
  sbatch -e $log -o $log -c 36 $batch $com
done
