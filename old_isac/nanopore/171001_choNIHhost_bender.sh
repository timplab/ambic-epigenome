#!/bin/bash
tag=171004_choNIHhost
rawroot=/dilithium/Data/Nanopore/oxford
wroot=/dilithium/Data/Nanopore/Analysis/
rawdir=$rawroot/$tag
wdir=$wroot/$tag
#mkdir $wdir
plotdir=/home/isac/Dropbox/Data/Nanopore/$tag/qc
#mkdir -p $plotdir


for tout in `ls $rawdir/*tout*`; do
  base=$(basename "$tout")
  base=${base%%.*}
  echo $tout
  echo $base
#  Rscript ~/Code/timp_nanopore/oxford/fq_check.R $tout $plotdir $base
#  ./plotCumLength.R $tout $plotdir $base
done
touts=`ls $rawdir/*tout*`
./combinedLengthPlot.R "$touts" $plotdir $tag


