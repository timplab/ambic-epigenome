#!/bin/bash
srcpath=/home/isac/Code/ilee/oxford
path=`readlink -f ./`
samp=$(basename "$path")
echo $samp
rawdir=/atium/Data/Nanopore/oxford/$samp
statpath=`ls ${rawdir}/*csv.gz`
outpath=${samp}_summary_stats.txt

${srcpath}/oxford_stats.R $statpath $outpath
