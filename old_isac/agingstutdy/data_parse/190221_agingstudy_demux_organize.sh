#!/bin/bash
root=/mnt
csv=$root/c/barcoding.csv
np=/home/ubuntu/Code/nanopolish/nanopolish

if [ "$1" == "index" ];then
  echo index
  for fq in $(find $root/*/reads/* -maxdepth 0 -type f -name "*fastq.gz"); do
    base=$(basename "$fq")
    rdir=$(dirname "$fq")
    base=${base%%.*}
    db=$(find $rdir/* -maxdepth 0 -type f -name "$base*readdb")
    if [ ! -z "$db" ];then
      continue
    fi
    pres="$pres $rdir/$base"
  done
  rdir={}
  fq={}.fastq.gz
  log={}.index.log
  com="$np index -d $rdir $fq &> $log"
  parallel $com ::: $pres
exit
fi

while IFS=$',' read -r -a line || [[ -n "$line" ]]
do
  runname=${line[0]}
  dir=$(find $root/*/demux/* -maxdepth 0 -type d -name "$runname")
  if [ -z $dir ];then
    continue
  fi
  barcode=${line[5]}
  echo $runname,$barcode
  grep $barcode $dir/fq_demux.log
  name=${line[1]}
  rdir=$(dirname "$(dirname $dir)")/reads
  # f5
  f5dir=$(find $dir/fast5 -maxdepth 1 -type d -name "$barcode")
  if [ -z $f5dir ];then
    continue
  fi
  newf5dir=$rdir/$name
  mv $f5dir $newf5dir
  # fq
  fq=$(find $dir/* -maxdepth 0 -type f -name "$barcode.fastq.gz")
  newfq=$rdir/$name.fastq.gz
  mv $fq $newfq

done <<< "$(awk 'NR>1' $csv)"
