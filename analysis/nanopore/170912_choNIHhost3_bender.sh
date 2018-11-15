#!/bin/bash
root=/atium/Data/Nanopore/Analysis
#srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
base=170912_choNIHhost3
echo $base
#tar=`ls ${tarpath}/${base}.raw.tgz` 
wpath=$root/$base
#mkdir -p $wpath
#tarpath=$wpath/tar


tout=`readlink -f ${wpath}/*csv.gz`

Rscript ~/Code/timp_nanopore/oxford/fq_check.R ${wpath}/*csv.gz ./ $base
