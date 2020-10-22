#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171024_cho
fqdir=${root}/fastq
bamdir=${root}/bam
[ -d $bamdir ] || mkdir $bamdir

align=/home/isac/Code/ilee/oxford/oxford_align.sh
#aligner=ngmlr
## try minimap2
aligner=minimap2
ref=/mithril/Data/NGS/Reference/cho/criGri1_plasmid2.fa

for fq in `find $fqdir -name "*fq.gz"`
do
  samp=$(basename "$fq")
  samp=${samp%.fq.gz}
  echo "aligning $samp"
  $align -i $fq -a $aligner -s $samp -r $ref -o $bamdir
done

