#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
logdir=$root/log/sv
[ -e $logdir ]||mkdir $logdir
bamdir=$root/pooled/bam
mbeddir=$root/pooled/methylation/methbyread
insertiondir=$root/insertion
[ -e $insertiondir ]||mkdir $insertiondir

for samp in Host StableGlut UnstableGlut;do
  label=${samp}D$day
  echo $label
  bam=`readlink -f $bamdir/$label.pooled.bam`
  insertion=$insertiondir/$label.insertion.sam
  insertnames=$insertiondir/$label.insertion.rnames.txt
  if [ ! -e $insertion ];then
    samtools view -h $bam \
      ambic_sigma_IgG_LC ambic_sigma_IgG_HC > $insertion
  fi
  samtools view $insertion | grep scaffold | cut -f1 > $insertnames
  #sort $insertnames | uniq -c
  insertmeth=$insertiondir/$label.insertion.methbyread.bed.gz
  mbed=`find $mbeddir -name "*$label*gz"`
  if [ ! -e $insertmeth ];then
    insertnum=`wc -l < $insertnames`
    if [ $insertnum -eq 0 ];then
      continue
    fi
    zcat $mbed |\
      grep -F -f $insertnames |\
      sort -k1,1 -k2,2n |\
      bgzip > $insertmeth
    tabix -p bed $insertmeth
  fi
done
