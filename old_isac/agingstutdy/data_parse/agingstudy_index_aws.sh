#!/bin/bash
TMPDIR=/mnt/c/tmp
export TMPDIR=$TMPDIR TMP=$TMPDIR TEMP=$TMPDIR

root=/mnt
np=/home/ubuntu/Code/nanopolish/nanopolish
#coms=$root/e/parallel/index_parallel
#[ ! -e $coms ]||rm $coms
for fq in $(find $root/*/reads -maxdepth 1 -name "*fastq.gz"); do
  pre=${fq%%.*}
  db=$fq.index.readdb
  if [ -e $db ]; then
    continue
  fi
  args="$args $pre"
# multi f5 + summary is not suppoorted
#  sum=$pre.summary.txt
#  if [ -e $sum ];then
#    arg="-s $sum"
#  fi
#  log=$pre.index.log
#  echo "$np index -d $pre -s $sum $fq &> $log" >> $coms
done
#parallel ::: < $coms
pre="{}"
fq=$pre.fastq.gz
log=$pre.index.log
com="$np index -d $pre $fq &> $log"
if [ "$1" == "dryrun" ];then
  com="echo {}"
fi
parallel $com ::: $args
