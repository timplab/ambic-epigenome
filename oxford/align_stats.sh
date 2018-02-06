#!/bin/bash

bampath=$1
bam=$(basename "$bampath")
prefix=${bam%%.*}

header=",aligned,unaligned,total,fraction"
line1=`samtools idxstats $bampath | awk 'BEGIN { OFS="," }{ al+=$3;un+=$4 } \
  END { all=al+un ; print "reads",al,un,all,al/all }'`
echo $header
echo $line1
line2=`awk 'BEGIN { OFS=",";FS="," }{ al+=$2 } \
  END { print al }'`
echo $line2




