#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
logdir=$root/log/sv
[ -e $logdir ]||mkdir $logdir
bamdir=$root/pooled/bam
svdir=$root/sv
[ -e $svdir ]||mkdir $svdir

snifflesarg="-s 2 -n 30"
sniffles=`readlink -f /home/isac/Code/Sniffles/bin/sniffles*/sniffles`

for samp in Host StableGlut UnstableGlut;do
  label=${samp}D$day
  bam=`readlink -f $bamdir/$label.pooled.bam`
  svout=$svdir/$label.sniffles.vcf
  echo $bam
  log=$logdir/$label.sniffles.log
  com="$sniffles $snifflesarg -m $bam \
    --tmp_file $svdir -v $svout >> $log"
  echo $com > $log
  eval $com
done
  

