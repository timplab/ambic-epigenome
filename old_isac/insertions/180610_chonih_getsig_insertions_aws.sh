#!/bin/bash
root=/shared/data/analysis
bamdir=$root/pooled/bam
insertiondir=$root/insertion
svdir=$root/sniffles
svouts=$(find $svdir -name "*NIH*.vcf")
plasmid=4_0cdhfr_vrc01wtg1m3_dgv
for svout in $svouts;do
  base=$(basename "$svout")
  label=${base%%.*}
  echo $label
  out=$insertiondir/$label.insertions.txt
#  grep 4_0 $svout | cut -f1,2,3,5,10 | grep -v 4_0 > $out
done

regname=scaf93
chrom=scaffold_93
bp=506559
region="$chrom:$upstream-$bp"
reg=/dilithium/Data/Nanopore/Analysis/171025_cho/insertion/picr49_insertion.bed
svout=$(find $svdir -name "*NIHIgG*vcf")

# get names of reads that support my favorite SV
rnames=$insertiondir/choNIHIgG.$regname.insertion.rnames.txt
#grep $chrom$'\t'$bp $svout | cut -d";" -f10 |\
#  sed 's/RNAMES=//' | tr "," "\n" > $rnames

# do this for all insetions
insertions=$insertiondir/choNIHIgG.insertions.txt
rnamesall=$insertiondir/choNIHIgG.insertions.rnames.txt
cat $insertions | while IFS=$'\t' read -r -a line
do
  chrom=${line[0]}
  bp=${line[1]}
#  grep $chrom$'\t'$bp $svout | cut -d";" -f10 |\
#    sed 's/RNAMES=//' | tr "," "\n" >> $rnamesall
done

# subset bam
bam=$(find $bamdir -name "*IgG*bam")
insertion=$insertiondir/choNIHIgG.insertion.sam
if [ ! -e $insertion ];then
  echo $bam
  samtools view $bam \
    $plasmid | grep -F -f $rnamesall  > $insertion
fi

