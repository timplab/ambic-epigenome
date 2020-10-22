#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171025_cho
mcalldir=$root/meth
outdir=$root/meth/readlevel

parse=/home/isac/Code/ilee/nanopolish/script/mtsv2bedGraph.py

for samp in host IgG;do
  for rep in 1 2 3;do
    label=choNIH${samp}$rep
    mbed=$outdir/$label.readlevel.meth.bed.gz
    tsv=$(find $mcalldir -name "$label*tsv.gz")
    com="gunzip -c $tsv | python $parse |\
      sort -T $outdir -k1,1 -k2,2n |\
      bgzip > $mbed ;\
      tabix -p bed $mbed"
    echo $com
    eval $com
  done
done
