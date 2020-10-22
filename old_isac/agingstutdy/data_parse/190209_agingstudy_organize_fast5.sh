#!/bin/bash
root=/mnt/sigma
root=/home/ambic
indir=$root/reads_unorg
outroot=$root/reads
[ -e $outroot ]||mkdir $outroot
scr=/home/ubuntu/Code/ilee/oxford/fast5_reorganize_renamedir.py
awk 'NR>1' $root/samples.csv | while IFS=$',' read -r -a line || [[ -n "$line" ]]
do
  args="-s ${line[2]}"
  base=${line[0]}
  dir=$(find $indir -maxdepth 1 -name "$base")
  if [ -z $dir ];then
    continue
  fi
  newname=${line[4]}
  if [[ $newname == *"PAD"* ]];then
    n=100000
  else
    n=10000
  fi
  newdir=$outroot/$newname
  slot=${line[3]}
  if [ ! -z "$slot" ];then
    args="$args -s $slot"
    exit
  fi
  com="python $scr -b $base $args -n $n -o $newdir $dir"
  echo $com
  eval $com
done
