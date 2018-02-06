#!/bin/bash

wdir=/dilithium/Data/NGS/Aligned/171025_choatac
bamdir=${wdir}/bam
qcdir=${wdir}/qc

## first do the idxstats
if [ 0 -eq 1 ];then
  for bam in `find $bamdir -name "*nodup.bam"`;do
    samp=$(basename "$bam")
    base=${samp%.*}
    echo $base
    samtools idxstats $bam > $bamdir/$base.idxstats.txt
  done
fi

## use idxstats to get out number of reads in the genome
if [ 0 -eq 1 ];then
  rnumout=$qcdir/align.readnum.txt
  echo "sample reads" | tr " " "\t" > $rnumout
  for stat in `find $bamdir -name "*nodup.idxstats.txt"`;do
    samp=$(basename "$stat")
    base=${samp%%.*}
    num=`awk '{ sum+=$3 }END{ print sum }' $stat`
    echo "$base $num" | tr " " "\t" >> $rnumout
  done
fi
