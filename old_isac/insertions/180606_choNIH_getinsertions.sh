#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
bamdir=$root/bam/pooled
insertiondir=$root/insertion

for samp in host IgG;do
  echo $samp
  bam=`find $bamdir -name "*$samp*bam"`
  insert=$insertiondir/$samp.insertion
  if [ ! -e $insert.sam ];then
    samtools view $bam |\
      grep 4_0cdhfr_vrc01wtg1m3_dgv > $insert.sam
  fi
    cut -f3 $insert.sam | sort | uniq -c | sort -k1,1n
done
