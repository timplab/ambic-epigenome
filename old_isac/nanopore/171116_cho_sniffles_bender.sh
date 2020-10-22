#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
########aligner=minimap2
bamdir=${root}/bam/pooled
svdir=${root}/sniffles
[ -e $svdir ]||mkdir $svdir
ref=/mithril/Data/NGS/Reference/cho/picr_IgG2/picr_IgG2.fa
srcdir=~/Code
sniffles=`readlink -f $srcdir/Sniffles*/bin/sniffles*/sniffles`
cfile=$root/sniffles.command.txt
cfile2=$root/sniffles.bedpe.command.txt
#rm $cfile

for bam in `find $bamdir -name "*.pooled.bam"`
do
  samp=$(basename "$bam")
  samp=${samp%%.*}
#  echo $samp
  bedpe=$svdir/$samp.sniffles.bedpe
  vcf=$svdir/$samp.sniffles.vcf
  bedpelog=$bedpe.log
  vcflog=$vcf.log
#  $sniffles -m $bam -b $bedpe -s 3 &> $bedpelog
#  $sniffles -m $bam -v $vcf -s 1 -n -1 &> $vcflog
  vcf=$svdir/$samp.sniffles.stringent.vcf
  vcflog=$vcf.log
#  echo "$sniffles -m $bam -v $vcf -s 5 -n -1 &> $vcflog" >> $cfile
  bedpe=$svdir/$samp.sniffles.stringent.bedpe
  bedpelog=$bedpe.log
#  echo "$sniffles -m $bam -b $bedpe -s 5 &> $bedpelog" >> $cfile2

done

parallel ::: < $cfile2
