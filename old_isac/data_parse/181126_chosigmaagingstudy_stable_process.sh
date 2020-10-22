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
  [ -e $logdir ]||mkdir -p $logdir
  for base in $samples; do
    fq=$(find $fqroot -type f -name "$base*fastq*" ! -name "*index*")
    bam=$(find $bamdir -name "$base.*sorted.bam")
    mout=$mdir/$base.meth.tsv
    log=$logdir/$base.mcall.log
    com="/shared/bin/nanopolish call-methylation -v -t 36 -r $fq -g $ref -b $bam > $mout"
    echo $com
    sbatch -e $log -o $log -c 36 -J "mcall" -t 100:0:0 $batch "$com"
  done
fi
if [ "$1" == "gzip" ];then
  for mout in $(find $mdir -name "*meth.tsv");do
    log=$mdir/$base.gzip.log
    if [ -e "$mout.gz" ];then
      echo exists,deleting
      rm $mout.gz
    fi
    echo $mout
    sbatch -e $log -o $log $batch "gzip $mout"
  done
fi

if [ "$1" == "fqgz" ];then
  for fq in $(find $fqroot -type f -name "*fastq");do
    echo $fq
    log=$logroot/$base.fqgzip.log
    sbatch -e $log -o $log $batch "gzip $fq"
  done
fi
