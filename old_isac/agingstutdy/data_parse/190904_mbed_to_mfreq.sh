#!/bin/bash
export LC_ALL=C
root=/mnt/d/Data/pooled_rep
bdir=$root/mbed
fdir=/mnt/e/Data/pooled_rep/mfreq
[ -e $fdir ]||mkdir $fdir -p

scr=/home/ubuntu/Code/nanopore-methylation-utilities/parseMethylbed.py

for mbed in $(find $bdir -name "*bed.gz"); do
  base=$(basename "${mbed%.bed.gz}")
  bases="$bases $base"
done

base="{}"
mbed=$bdir/$base.bed.gz
mfreq=$fdir/$base.freq.txt.gz
com="gunzip -c $mbed | $scr frequency |\
  gzip > $mfreq"

parallel $com ::: $bases

