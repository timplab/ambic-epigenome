#!/bin/bash

script=./script/atacseqPipeline.py
wdir=/dilithium/Data/NGS/Aligned/180223_choOmniAtac
btidx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd
kdir=${wdir}/180223_choOmniAtac.metadata.csv
plotdir="/home/isac/Dropbox/Data/Atac-seq/180223_choOmniAtac/plot"

#$script -k $kdir -o $wdir -x $btidx -g mm -m bamprocess

beddir=$wdir/bed
if [ "$1" == "simulate" ];then
  python ~/Code/ilee/util/simulate_random_coverage.py /mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa.fai -d 0.01 | sort -k1,1 -k2,2n | gzip > $beddir/simulate.bed.gz
fi

tss=/mithril/Data/NGS/Reference/cho/chok1/annotation/Cricetulus_griseus_chok1gshd.TSS.bed
tssdir=$wdir/tss

if [ "$1" == "getdist" ];then
  for bed in `find $beddir -name "*bed.gz"`;do
    base=$(basename "$bed")
    base=${base%%.*}
    echo $base
    gunzip -c $bed |\
      bedtools closest -D b -b $tss -a stdin |\
      awk '{ if(($NF>-2000)&&($NF<2000)) print }'> $tssdir/$base.tss.bed
    # parse
    awk '{ OFS="\t" }{ print $NF,$10 }' $tssdir/$base.tss.bed > $tssdir/$base.tssdist.txt
  done
fi

if [ "$1" == "plot" ];then
  for bed in `find $beddir -name "*bed.gz"`;do
    base=$(basename "$bed")
    base=${base%%.*}
    echo $base
    Rscript ~/Code/ambic-epigenome/plot/atacSeq_plots.R covByDistance -i $tssdir/$base.tssdist.txt -o $plotdir/$base.tssdist.pdf
  done
fi
