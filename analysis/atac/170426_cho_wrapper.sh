#!/bin/bash
qccodedir=/home/isac/Code/ilee/qc
ref=/mithril/Data/NGS/Reference/cho/picr.corrected.fa
btidx=/mithril/Data/NGS/Reference/cho/cho_picr
rawdir=/atium/Data/NGS/Raw/170426_ILSW_CHO
bamdir=/atium/Data/NGS/Aligned/170426_cho/bam
qcdir=/atium/Data/NGS/Aligned/170426_cho/qc
base=choIgGATAC

for i in 1 2 3; do
  samp=${base}${i}
  echo $samp
  r1=`ls ${rawdir}/${samp}*R1*fastq.gz`
  r2=`ls ${rawdir}/${samp}*R2*fastq.gz`
  sam=${bamdir}/${samp}.sam
  bam=${bamdir}/${samp}.bam
  sorted=${bamdir}/${samp}.sorted.bam
  
  #fastq QC
  #${qccodedir}/fastq_qc.sh $r1 $qcdir
  #${qccodedir}/fastq_qc.sh $r2 $qcdir


  ## align
  #bowtie2 -X 2000 -p 10 --time \
  #  -x $btidx -1 $r1 -2 $r2 \
  #  -S $sam &> ${bamdir}/${samp}_bowtie2.log

  ## samtools
  #samtools view -bS $sam > $bam
  #samtools sort $bam -o $sorted 
  #samtools index $sorted
  #rm $sam $bam

  ##Picard

  ## bam qc
  #${qccodedir}/bam_qc.sh $sorted $qcdir

  
done

#cat ${qcdir}/*fastqQC.csv > ${qcdir}/${base}_fastqQCall.csv
#cat ${qcdir}/*bamQC.csv > ${qcdir}/${base}_bamQCall.csv





