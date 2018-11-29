#!/bin/bash
root=/shared/Data
demuxdir=$root/demux
fqroot=$root/fastq
logdir=$root/log/consoldiate
[ -e $logdir ]||mkdir $logdir
batch=/shared/Code/ilee/slurm/batchcommand.scr

if [ "$1" == "consolidate" ];then
  bcsv=$(find $root -maxdepth 1 -name "*barcodes.csv")
  awk 'NR>1' $bcsv | while IFS=$',' read -r -a line
  do
    barcode=${line[0]}
    lab=${line[9]}
    run=${line[1]}
    if [ "${line[8]}" != "rapid" ];then
      continue
    fi
    if [ -z "$barcode" ];then
      continue
    fi
    if [ $barcode -lt 10 ];then
      barcode=0$barcode
    fi
    echo $lab
    fqname=BC$barcode.fastq
    fqdir=$(find $demuxdir -type d -name "$run")
    outdir=${fqdir#$demuxdir/}
    outdir=$fqroot/$(dirname "$outdir")
    outfq=$outdir/$lab.fastq
    fqs=$(find $fqdir -name $fqname)
    com="cat $fqs > $outfq"
    log=$logdir/$lab.consolidate.log
    sbatch -o $log -e $log -c 1 -J con-$lab $batch "$com"
  done
fi

