#!/bin/bash
awsroot="s3://timp.ambic/choSigmaAgingStudy/nanopore"
dir=/mnt/d/Data/allbam
bdir=/mnt/d/Data/bed
qc=$bdir/bedqc.csv
echo "sample,numreads,numbases,n50,meanlen,q1,median,q3,max" > $qc
day=0
for samp in Host StableGln StableNogln UnstableGln UnstableNogln; do
  base=CHOZN${samp}Day$day
  if [ "$1" == "mergebam" ];then
    for rep in 1 2 3; do
      pre=${base}_$rep
      echo $pre
      bams=$(find $dir -name "$pre*sorted.bam")
      samtools merge -@ 36 $dir/$pre.bam $bams
      samtools index $dir/$pre.bam
      rm $bams
      aws s3 sync $dir/ $awsroot/pooled_rep/ --exclude "*" --include "$pre.bam" --include "$pre.bam.bai"
    done
    echo $base
    bams=$(find $dir -name "$base*bam")
    samtools merge -@ 36 $dir/$base.pooled.bam $bams
    samtools index $dir/$base.pooled.bam
    rm $bams
    aws s3 sync $dir/ $awsroot/pooled/ --exclude "*" --include "$base.pooled.bam" --include "$base.pooled.bam.bai"
  fi
  if [ "$1" == "bamtobed" ];then
    bedtools bamtobed -i $dir/$base.pooled.bam  > $bdir/$base.pooled.bed
  fi
  if [ "$1" == "bedqc" ];then
    python /home/ubuntu/Code/ilee/qc/bedQC.py $bdir/$base.pooled.bed >> $qc
  fi


done

