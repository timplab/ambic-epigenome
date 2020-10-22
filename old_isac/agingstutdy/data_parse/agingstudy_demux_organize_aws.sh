#!/bin/bash
root=/mnt/e/
sroot=/mnt/e/reads_split
f5droot=/mnt/e/reads_demux
outroot=$root/reads
demuxdir=$root/demux
csv=$demuxdir/samples.csv

while IFS=$',' read -r -a line || [[ -n "$line" ]]
do
  runname=${line[0]}
  barcode=${line[5]}
  sample=${line[1]}
  outdir=$outroot/$sample
  # first fastq
  bcdir=$demuxdir/$runname
  if [ "$1" == "fq" ];then
    fqdir=$bcdir/$barcode
    args="$args $fqdir $sample"
  fi
  if [ -e $outdir ];then
    continue
  fi
  echo $sample,$barcode
  if [ "$1" == "f5" ];then
    sdir=$sroot/$runname
    f5dir=$f5droot/$runname
    if [ ! -e $sdir ];then
      dir=$(find $root/reads_unorg/* -maxdepth 0 -type d -name "$runname")
      multi_to_single_fast5 --recursive -t 36 \
        -i $dir -s $sdir
      rm -r $dir/*
    fi
    if [ ! -e $f5dir ];then
      bsum=$bcdir/barcoding_summary.txt
      python /home/ubuntu/Code/ilee/oxford/fast5_split_barcodes.py \
        -o $f5dir -b $bsum -v -i $sdir
    fi
    [ -e $outdir ]||mkdir $outdir
    single_to_multi_fast5 --recursive -t 36 -n 100000 \
      -i $f5dir/$barcode -s $outdir -f $sample 
  fi
done <<< "$(cat $csv)"
if [ "$1" == "fq" ];then
  fqs={1}/*fastq
  outfq=$outroot/{2}.fastq.gz
  com="cat $fqs | gzip > $outfq"
#  com="echo {1},{2}"
  parallel -j 36 -N 2 $com ::: $args
fi
