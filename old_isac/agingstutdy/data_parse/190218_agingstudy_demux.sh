#!/bin/bash
TMPDIR=/mnt/c/tmp
TMP=$TMPDIR
TEMP=$TMPDIR
export TMP TEMP TMPDIR
bcodedir=/mnt/c/barcode
samples=/mnt/c/samples.csv
samps=CHOZNUnstableNoglnDay6090_3_PAD09846
scr=/home/ubuntu/Code/ilee/oxford/fast5_split_barcodes.py

while IFS=$',' read -r -a line || [[ -n "$line" ]]
do
  base=${line[4]}
  barcode=${line[7]}
  if [ "$barcode" != "Y" ];then
    continue
  fi
  bc=$(find $bcodedir/$base -name "barcoding_summary.txt")
  if [ -z $bc ]; then
    continue
  fi
  fq=$(find /mnt/*/reads/* -maxdepth 0 -type f -name "$base.fastq.gz")
  if [ ! -z $fq ];then
    root=$(dirname "$(dirname $fq)")
    outdir=$root/demux/$base
    args="$args -b $bc -f $fq -o $outdir $outdir/fq_demux.log"
  fi
  db=$(find /mnt/*/reads/* -maxdepth 0 -type f -name "$base.*readdb")
  if [ ! -z $db ];then
    root=$(dirname "$(dirname $db)")
    args="$args -b $bc -r $db -o $outdir $outdir/f5_demux.log"
  fi
  if [ ! -z $root ]&&[ "$1" == "clean" ];then
    echo cleaning $base
    mv $root/reads/$base* $root/reads_unused/
  fi
done <<< "$(awk 'NR>1' $samples)"
if [ "$1" == "clean" ];then
  exit
fi
com="python -u $scr -v {1} {2} {3} {4} {5} {6} &> {7}"
#com="echo {1} {2} {3} {4} {5} {6} {7}"
parallel -N 7 $com ::: $args


