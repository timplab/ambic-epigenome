#!/bin/bash
qccodedir=/home/isac/Code/ilee/qc
root=/atium/Data/NGS/Aligned/170718_ambicatac
bamdir=${PWD}/bam
qcdir=${PWD}/qc
regionbed=${PWD}/insertions.bed

for samp in `ls ${bamdir}/*bam`; do
  echo $samp
  base=${samp##*/}
  base=${base%%.*}
  echo $base
  ##first getting chr stats
  idxstat=$qcdir/${base}_idxstats.txt
  #samtools idxstats $samp > $idxstat
  ##MT percentage
  echo "mapped,mito"
  echo `grep "NW" $idxstat | awk '{ sum += $3 } END { print sum }'`,`grep "NC" $idxstat | awk '{ print $3 }'`
  ## samtools to get coverage in interested regions
  samtools bedcov $regionbed $samp
  #samtools depth -a -b $regionbed $samp

done
