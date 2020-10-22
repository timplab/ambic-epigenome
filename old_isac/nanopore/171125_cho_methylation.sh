#!/bin/bash
wroot=/dilithium/Data/Nanopore/Analysis/171025_cho
methdir=$wroot/meth
runpath=$wroot/methfreq.txt
rm $runpath

for f in `find $methdir -maxdepth 1 -name "*meth.tsv.gz"`;do
  base=$(basename "$f")
  samp=${base%.meth.tsv.gz}
  echo $samp
  out=$methdir/$samp.methfreq.tsv
  echo "gunzip -c $f | \
    python /home/isac/Code/nanopolish/scripts/calculate_methylation_frequency.py \
    -c 2.5 > $out" >> $runpath
done
parallel ::: < $runpath
