#!/bin/bash
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/barcodingtest
f5dir=$root/fast5
bcalldir=$root/basecalled
[ -e $bcalldir ]||mkdir $bcalldir
dbindir=$root/deepbin
[ -e $dbindir ]||mkdir $dbindir

if [ "$1" == "bcall" ];then
  read_fast5_basecaller.py -k SQK-LSK109 -f FLO-PRO001 --barcoding \
    -r --disable_filtering -t 8 -s $bcalldir -i $f5dir
fi
if [ "$1" == "deepbinner" ];then
  com="deepbinner classify --verbose --native $f5dir/10"
  echo $com
  eval $com
fi
if [ "$1" == "summarize" ];then
  summary=$bcalldir/sequencing_summary.txt
  bsum=$bcalldir/barcoding_summary.txt
  python ./summarize_barcode.py $summary > $bsum
fi

