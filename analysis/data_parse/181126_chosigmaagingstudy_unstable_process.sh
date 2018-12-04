#!/bin/bash
root=/shared/Data
fqroot=$root/fastq
logroot=$root/log
bamdir=$root/bam
mdir=$root/mcall
[ -e $mdir ]||mkdir $mdir
batch=/shared/Code/ilee/slurm/batchcommand.scr
ref=/shared/Data/Reference/Cricetulus_griseus_chok1gshd.sigmaIgG.fa

npwrapper=/shared/Code/ilee/nanopolish/slurm/nanopolishWrapper.sh
np=/shared/bin/nanopolish

samples=$(awk '$2=="x"{ print $1 }' $root/samples.txt | tr "\n" " ")
if [ "$1" == "mcall" ];then
  logdir=$logroot/mcall
  [ -e $logdir ]||mkdir $logdir
  for base in $samples; do
    fq=$(find $fqroot -type f -name "$base*fastq*" ! -name "*index*")
    bam=$(find $bamdir -name "$base.*sorted.bam")
    mout=$mdir/$base.meth.tsv
    log=$logdir/$base.mcall.log
    com="/shared/bin/nanopolish call-methylation -v -t 36 -r $fq -g $ref -b $bam > $mout"
    echo $com
    sbatch -e $log -o $log -c 36 -t 100:0:0 $batch "$com"
  done
fi

if [ "$1" == "gzip" ];then
  for mtsv in $(find $mdir -name "*tsv"); do
    log=$mtsv.gzip.log
    com="gzip $mtsv"
    sbatch -e $log -o $log $batch "$com"
  done
fi

if [ "$1" == "fqgz" ];then
  for fq in $(find $fqroot -type f -name "*fastq"); do
    echo $fq
    log=$fq.gzip.log
    com="gzip $fq"
    sbatch -e $log -o $log $batch "$com"
  done
fi
