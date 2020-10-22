#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
bamdir=$root/bam/pooled
insertiondir=$root/insertion/picr
svout=/dilithium/Data/Nanopore/Analysis/171025_cho/sniffles/picr/choNIHIgG.sniffles.vcf
plasmid=4_0cdhfr_vrc01wtg1m3_dgv
grep 4_0 $svout | cut -f1,2,3,5,10 | grep -v 4_0

regname=picr49
chrom=picr_49
upstream=6218858
bp=6220858
region="$chrom:$upstream-$bp"
reg=/dilithium/Data/Nanopore/Analysis/171025_cho/insertion/picr49_insertion.bed

# get names of reads that support my favorite SV
rnames=$insertiondir/$regname.insertion.rnames.txt
grep $chrom$'\t'$bp $svout | cut -d";" -f10 |\
  sed 's/RNAMES=//' | tr "," "\n" > $rnames

bamdir=$root/bam/pooled
bam=$(find $bamdir -name "*IgG*bam")
insertion=$insertiondir/choNIHIgG.$regname.insertion.sam
if [ ! -e $insertion ];then
  echo $bam
  samtools view $bam \
    $plasmid | grep -F -f $rnames  > $insertion
fi

mbeddir=$root/meth/readlevel
# subset methylation by read name
for samp in IgG;do
  label=choNIH$samp
  mbed=$(find $mbeddir -name "$label*.bed.gz")
  insertmeth=$insertiondir/$label.$regname.insertion.meth.bed.gz
  if [ ! -e $insertmeth ];then
    zcat $mbed |\
      grep -F -f $rnames |\
      sort -k1,1 -k2,2n |\
      bgzip > $insertmeth
    tabix -p bed $insertmeth
  fi
done



