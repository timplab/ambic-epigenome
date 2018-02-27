#!/bin/bash

script=./script/atacseqPipeline.py
wdir=/dilithium/Data/NGS/Aligned/171025_choatac
beddir=$wdir/bed
tssdir=$beddir/tss
btidx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd
kdir=${wdir}/choatac.metadata.csv
plotdir=~/Dropbox/Data/Atac-seq/171025_choatac/plot

if [ "$1" == "trim" ];then
  $script -k $kdir -o $wdir -x $btidx -g mm -m trim
fi
if [ "$1" == "align" ];then
  $script -k $kdir -o $wdir -x $btidx -g mm -m align
fi
if [ "$1" == "process" ];then
  $script -k $kdir -o $wdir -x $btidx -g mm -m bamprocess
fi

if [ "$1" == "getdist" ] || [ "$1" == "plotdist" ];then
  echo "$1"
  ./script/atacCovByDistance.sh $1 -i $beddir -a /mithril/Data/NGS/Reference/cho/chok1/annotation/Cricetulus_griseus_chok1gshd.TSS.bed -t tss
fi
if [ "$1" == "plotmove" ];then
  cp $tssdir/*pdf $plotdir
fi
