#!/bin/bash
rawroot=/atium/Data/Nanopore/oxford
root=/home/isac/Dropbox/Data/Nanopore
base1=170922_choNIHhost1
base2=170922_choNIHhost2


for base in $base1 $base2; do
  wpath=$root/$base
  mkdir -p $wpath
  tout=`readlink -f $rawroot/$base/$base.tout.csv.gz`
  Rscript ~/Code/timp_nanopore/oxford/fq_check.R $tout $wpath $base
done

