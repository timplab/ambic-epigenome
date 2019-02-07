#!/bin/bash
srcdir=$(dirname $(readlink -f $0))
fa="$1"
out="$2"
binsize="$3"

if [ ! -e $fa.fai ];then
  samtools faidx $fa
fi

python $srcdir/../util/bin_genome.py -f $fa.fai -b $binsize > $out

