#!/bin/bash
# usage: ./bcall_yieldcheck.sh /path/to/summary

dir=$1
sums=`find $dir -name "*summary.csv.gz"`

gunzip -c $sums | awk '{ FS="," }{ sum+=$13 }END{ print sum/1000000000,NR,sum/NR }'
