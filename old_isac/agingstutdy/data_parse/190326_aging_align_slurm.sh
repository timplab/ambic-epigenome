#!/bin/bash
root=/shared/Data
fqdir=$root/fastq
bamdir=$root/bam
scr=/home/ubuntu/Code/ilee/oxford/slurm/oxford_align.scr
ref=/shared/Data/ref/Cricetulus_griseus_chok1gshd.sigmaIgG.fa
if [ "$1" == "upload" ];then
  while :
  do
    for log in $(find $bamdir -name "*log"); do
      lastline=$(tail -n1 $log)
      if [[ $lastline == "Finished"* ]];then
        base=$(basename "$log")
        base=${base%%.*}
        aws s3 sync $bamdir s3://timp.ambic/choSigmaAgingStudy/nanopore/bam/ --exclude "*" --include "$base*"
      fi
    done
  done
  exit
fi
for fq in $(find $fqdir -name "*fastq.gz"); do
  base=$(basename $fq)
  base=${base%%.*}
  bam=$bamdir/$base.sorted.bam
  log=$bamdir/$base.align.log
  if [ -e $log ];then
    echo "$base already aligning"
    continue
  fi
  echo $base

  sbatch -e $log -o $log -c 72 $scr \
    -r $ref -i $fq -a minimap2 -o $bamdir -e aws --samarg "-q 20"
done

