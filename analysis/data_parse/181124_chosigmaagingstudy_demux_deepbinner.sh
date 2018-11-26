#!/bin/bash
root=/shared/Data
tardir=$root/tar
logdir=$root/log/demux
[ -e $logdir ]||mkdir -p $logdir
fqroot=$root/fastq
f5root=$root/fast5
demuxroot=$root/demux_deepbinner
batch=/shared/Code/ilee/slurm/batchcommand.scr
fqs=$(find $fqroot -name "*f*q.gz")

n=1
tensorargs="--intra_op_parallelism_threads $n --omp_num_threads $n"
for fq in $fqs; do
  dir=$(dirname "$fq")
  dir=${dir#$fqroot/}
  base=$(basename "$fq")
  base=${base%%.*}
  if [ "$dir" == "D30" ];then
    continue
  fi
  tar=$(find $tardir/$dir -name "$base.*tgz")
  if [ $tar ];then
    # tar still there means untar hasn't been done/complete
    continue
  fi
  echo $base
  outdir=$demuxroot/$dir/$base
  [ -e $outdir ]||mkdir -p $outdir
  runtype="--native"
  if [ "$base" == "181106_choSigmaStableGlutD306090rep1" ];then
    # rapid
    runtype="--rapid"
  fi
  f5dir=$f5root/$dir/$base
  class=$outdir/$base.classifications
  com="/shared/bin/deepbinner classify $tensorargs $runtype $f5dir > $class"
  echo $com
  eval $com
done
  
