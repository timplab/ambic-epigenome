#!/bin/bash
export PATH="/shared/bin:$PATH"
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

if [ "$1" == "mcall" ];then
  logdir=$logroot/mcall
  [ -e $logdir ]||mkdir $logdir
  for fq in $(find $fqroot -type f -name "*fastq*" ! -name "*index*");do
    base=$(basename "$fq")
    base=${base%%.*}
    bam=$(find $bamdir -name "$base.*sorted.bam")
    mout=$mdir/$base.meth.tsv
    log=$logdir/$base.mcall.log
    com="nanopolish call-methylation -v -t 36 -r $fq -g $ref -b $bam > $mout"
    echo $com
    sbatch -e $log -o $log -c 36 -t 100:0:0 $batch "$com"
  done
fi
if [ "$1" == "gzip" ];then
  for mcall in $(find $mdir -name "*tsv");do
    echo $mcall
    log=$mcall.gzip.log
    sbatch -e $log -o $log $batch "gzip $mcall"
  done
fi
if [ "$1" == "fqgzip" ];then
  logdir=$logroot/fqgzip
  [ -e $logdir ]||mkdir $logdir
  for fq in $(find $fqroot -type f -name "*fastq");do
    echo $fq
    base=$(basename "$fq")
    log=$logdir/$base.gzip.log
    sbatch -e $log -o $log $batch "gzip $fq"
  done
fi
