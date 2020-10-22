#!/bin/bash
qccodedir=/home/isac/Code/ilee/qc
ref=/mithril/Data/NGS/Reference/cho/choMT.fasta
rawdir=/atium/Data/NGS/Raw/170426_ILSW_CHO
outdir=/atium/Data/NGS/Aligned/170426_cho
blastdir=${outdir}/blast
fastadir=${outdir}/fasta
base=choIgGATAC

## make database
#makeblastdb -in $ref -parse_seqids -dbtype nucl

## make fasta from fastq
fasta=${fastadir}/170426_choATAC_aligned.fasta
tag=170426_choATAC_aligned

#rm $fasta
#
#for i in 1 2 3; do
#  samp=${base}${i}
#  echo $samp
#  r1=`ls ${rawdir}/${samp}*R1*fastq.gz`
#  r2=`ls ${rawdir}/${samp}*R2*fastq.gz`
#  seqtk seq -A $r1 >> $fasta
#  seqtk seq -A $r2 >> $fasta
#  
#done

blastn -db $ref -query $fasta -outfmt 6 -out ${blastdir}/${tag}_blastn.txt






