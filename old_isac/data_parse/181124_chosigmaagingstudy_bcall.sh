#!/bin/bash
root=/shared/Data
fqroot=$root/fastq
f5root=$root/fast5
bcallroot=$root/basecall

samples="181106_choSigmaStableGlutD306090rep1 180914_choSigmaHostD30rep2_min "
for sample in $samples;do
  echo $sample
  f5dir=$(find $f5root -maxdepth 2 -name "$sample")
  bcalldir=$bcallroot/$sample
  if [ "$1" == "bcall" ];then
    args="-f FLO-MIN106"
    if [ "$sample" == "181106_choSigmaStableGlutD306090rep1" ];then
      # rapid
      args="$args -k SQK-RAD004 --barcoding"
    else
      args="$args -k SQK-LSK109"
    fi
    /shared/Code/ilee/oxford/slurm/bcallWrapper.sh basecall \
      -o $bcalldir -r $f5dir -b $sample -e aws -a "$args"
  fi
  if [ "$1" == "org" ];then
    fqs=$(find $bcalldir -type f -name "*fastq")
    pre=$fqroot/${f5dir#$f5root/}
    fq=$pre.fastq.gz
    sum=$pre.summary.txt
    sums=$(find $bcalldir -type f -name "sequencing_summary.txt")
    cat $fqs | gzip > $fq
    awk 'FNR != 1 || NR == 1{ print }' $sums > $sum
  fi
  
done
  
