#!/bin/bash
wroot=/dilithium/Data/Nanopore/Analysis/171025_cho
ins=$wroot/insertion
meth=$ins/meth
rldir=$meth/readlevel

for f in `find $meth -maxdepth 1 -name "*meth.tsv"`;do
  base=$(basename "$f")
  samp=${base%.meth.tsv}
  echo $samp
  out=$rldir/$samp.readlevel.meth.tsv
#  ~/Code/ilee/nanopolish/readlevel_methylation.py -i $f > $out
  grep "picr" $out
done
