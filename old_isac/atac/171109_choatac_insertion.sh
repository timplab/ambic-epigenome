#!/bin/bash
root=/dilithium/Data/NGS/Aligned/171025_choatac
beddir=$root/bed
insertdir=$root/insertion
insertbed=$insertdir/choatac.insertion.bed

for samp in host IgG
do
  beds=`ls $beddir/*$samp*pool.bed.gz`
  echo $samp
  gunzip -c $beds | bedtools coverage -a $insertbed -b stdin > ${insertdir}/$samp.insert.10kb.bed
done
