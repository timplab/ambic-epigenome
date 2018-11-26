#!/bin/bash
root=/shared/Data
logdir=$root/log/demux
[ -e $logdir ]||mkdir -p $logdir
chunkroot=$root/fastq_chunks
demuxroot=$root/demux
batch=/shared/Code/ilee/slurm/batchcommand.scr
bases="181112_choSigmaStableGlutD306090rep1 181117_choSigmaStableNoGlutD0rep123 181117_choSigmaStableNoGlutD306090rep1"

for base in $bases; do
  dir=$(find $chunkroot -type d -name "$base")
  dir=$(dirname $dir)
  dir=${dir#$chunkroot/}
  chunkdir=$chunkroot/$dir
  pre=$chunkdir/$base/$base.fastq
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
done

