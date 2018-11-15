#!/bin/bash
tag=choIgGNIH
rawdir=/dilithium/Data/Nanopore/oxford/171025_cho
wdir=/dilithium/Data/Nanopore/Analysis/171025_cho/qc

dropdir=/home/isac/Dropbox/Data/ambic/nanopore/plots


for fh in `ls $rawdir/*$tag*csv*`; do
  base=$(basename "$fh")
  base=${base%%.*}
  if [ ! -e $wdir/${base}_fq.pdf ];then
    echo "$base"
    Rscript ~/Code/timp_nanopore/oxford/fq_check.R $fh $wdir $base
  fi
#  ./plotCumLength.R $samp $plotdir $base
done
#touts=`ls $rawdir/*tout*`
#./combinedLengthPlot.R "$touts" $plotdir $tag


