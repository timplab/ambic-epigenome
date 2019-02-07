#!/bin/bash
qccodedir=/home/isac/Code/ilee/qc
srcdir=/home/isac/Code/timp_genetics/atac
ref=/mithril/Data/NGS/Reference/cho/picr.corrected.fa
btidx=/mithril/Data/NGS/Reference/cho/GCF_000223135.1_CriGri_1.0_genomic
rawdir=/atium/Data/NGS/Raw/170824_BAILRW
outdir=/atium/Data/NGS/Aligned/170824_chok1ATAC
bamdir=$outdir/bam
qcdir=$outdir/qc
trimdir=$outdir/trimfq
base=chok1atac

for i in `ls ${rawdir}/${base}*_R1_*fastq.gz`; do
  samp=$(basename "$i")
  samp=${samp%%_*}
  echo $samp
  r1=`ls ${rawdir}/${samp}*R1*fastq.gz`
  r2=`ls ${rawdir}/${samp}*R2*fastq.gz`
  sam=${bamdir}/${samp}.sam
  bam=${bamdir}/${samp}.bam
  sorted=${bamdir}/${samp}.sorted.bam
  
  #fastq QC
  #${qccodedir}/fastq_qc.sh $r1 $qcdir
  #${qccodedir}/fastq_qc.sh $r2 $qcdir
  #fastqc -o $qcdir $r1 $r2

  ##trim
  #~/Code/trim_galore_zip/trim_galore -q 28\
  #   -o ${trimdir} --paired $r1 $r2 &> $trimdir/${samp}_trim.log
  trim1=`ls ${trimdir}/${samp}*val_1.fq.gz`
  trim2=`ls ${trimdir}/${samp}*val_2.fq.gz`

  ## align
  #$srcdir/atacAlign.sh -x $btidx -s $samp -t 4 -o $bamdir $trim1 $trim2

  ## bam qc
#  ${qccodedir}/bam_qc.py $sorted  > ${qcdir}/${samp}_bamQC.csv
  
done

awk '{OFS=","}{print FILENAME, $0}' ${qcdir}/*bamQC.csv \
  > ${qcdir}/chok1atac_all_bamQc.csv

#cat ${qcdir}/*fastqQC.csv > ${qcdir}/${base}_fastqQCall.csv
#cat ${qcdir}/*bamQC.csv > ${qcdir}/${base}_bamQCall.csv





