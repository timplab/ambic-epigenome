#!/bin/bash
root=/shared/Data
demuxdir=$root/demux
f5root=$root/fast5
fqroot=$root/fastq
logroot=$root/log
bamdir=$root/bam
batch=/shared/Code/ilee/slurm/batchcommand.scr

npwrapper=/shared/Code/ilee/nanopolish/slurm/nanopolishWrapper.sh
np=/shared/bin/nanopolish

if [ "$1" == "index" ];then
  logdir=$logroot/index
  [ -e $logdir ]||mkdir $logdir
  bcsv=$(find $root -maxdepth 1 -name "*barcodes.csv")
  awk 'NR>1' $bcsv | while IFS=$',' read -r -a line
  do
    lab=${line[9]}
    run=${line[1]}
    fq=$(find $fqroot -name "$lab.fastq*")
    sum=$(find $fqroot -name "$run.summary.txt")
    f5dir=$(find $f5root -maxdepth 2 -type d -name "$run")
    com="$np index -v -d $f5dir -s $sum $fq"
    log=$logdir/$lab.index.log
    echo $com
    sbatch -e $log -o $log -J "index-$lab" $batch "$com"
  done
fi
if [ "$1" == "align" ];then
  mmi=/shared/Data/Reference/Cricetulus_griseus_chok1gshd.sigmaIgG.mmi
  [ -e $bamdir ]||mkdir $bamdir
  logdir=$logroot/align
  [ -e $logdir ]||mkdir $logdir
  fqs=$(find $fqroot -type f -name "*fastq*" ! -name "*index*")
  for fq in $fqs; do
    base=$(basename "$fq")
    base=${base%%.*}
    pre=$bamdir/$base.chok1gshd_sigmaIgG
    outbam=$pre.bam
    outsam=$pre.sam
    sortbam=$pre.sorted.bam
    com="/shared/bin/minimap2 -ax map-ont -t 36 -L $mmi $fq > $outsam"
    log=$logdir/$base.align.log
    j=align
    if [ -e $sortbam ];then
      com="/shared/bin/samtools index $sortbam"
      log=$logdir/$base.index.log
      j=index
    elif [ -e $outbam ];then
      com="/shared/bin/samtools sort -T $pre.sorting -o $sortbam $outbam"
      log=$logdir/$base.sort.log
      j=sort
    elif [ -e $outsam ];then
      com="/shared/bin/samtools -@ 35 view -q 20 -b $outsam > $outbam"
      log=$logdir/$base.bam.log
      j=bam
      tail -n1 $logdir/$base.align.log
    fi
    echo $com
    sbatch -t 24:0:0 -c 36 -e $log -o $log -J $j $batch "$com"
  done
fi

