#!/bin/bash
root=/shared/Data
logdir=$root/log/demux
[ -e $logdir ]||mkdir -p $logdir
fqroot=$root/fastq
chunkroot=$root/fastq_chunks
demuxroot=$root/demux
batch=/shared/Code/ilee/slurm/batchcommand.scr
fqs=$(find $fqroot -name "*f*q.gz")

for fq in $fqs; do
  dir=$(dirname "$fq")
  dir=${dir#$fqroot/}
  base=$(basename "$fq")
  base=${base%%.*}
  if [ "$dir" == "D30" ];then
    continue
  fi
  size=$(du $fq | cut -f1)
  if [[ $size -gt 30000000 ]];then
    echo "split $base into 500e5 reads"
    chunkdir=$chunkroot/$dir/$base
    [ -e $chunkdir ]||mkdir -p $chunkdir
    n=2000000
    pre=$chunkdir/$base.fastq
    gunzip -c $fq | split -l $n - $pre
    for sub in $(find $chunkdir -name "$base.fastq*");do
      ext=${sub#$pre}
      echo $ext
      outdir=$demuxroot/$dir/$base/$ext
      log=$logdir/$base.$ext.demux.log
      [ -e $outdir ]||mkdir -p $outdir
      com="/shared/bin/porechop --threads 36 -i $sub -b $outdir"
      echo $com
      sbatch  -o $log -e $log -c 36 -J demux-$base $batch "$com"
    done
    continue
  fi
  outdir=$demuxroot/$dir/$base
  log=$logdir/$base.demux.log
  [ -e $outdir ]||mkdir -p $outdir
  com="/shared/bin/porechop --threads 36 -i $fq -b $outdir"
  echo $com
  sbatch  -o $log -e $log -c 36 -J demux-$base $batch "$com"
done
  
