#!/bin/bash
tag=171021_choNIHhost
rawroot=/dilithium/Data/Nanopore/oxford
wroot=/dilithium/Data/Nanopore/Analysis/
rawdir=$rawroot/$tag
wdir=$wroot/$tag
mkdir $wdir
plotdir=/home/isac/Dropbox/Data/Nanopore/$tag/qc
mkdir -p $plotdir


for samp in `ls $rawdir/*summary*`; do
  base=$(basename "$samp")
  base=${base%%.*}
  echo $samp
  echo $base
  Rscript ~/Code/timp_nanopore/oxford/fq_check.R $samp $plotdir $base
#  ./plotCumLength.R $samp $plotdir $base
done
#touts=`ls $rawdir/*tout*`
#./combinedLengthPlot.R "$touts" $plotdir $tag


