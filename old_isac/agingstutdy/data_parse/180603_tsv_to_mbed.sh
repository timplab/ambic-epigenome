#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
pooldir=$root/pooled/methylation/methbyread
[ -e $pooldir ]||mkdir -p $pooldir
mcalldir=$root/mcall-cpg

parse=/home/isac/Code/ilee/nanopolish/script/mtsv2bedGraph.py

for samp in Host StableGlut UnstableGlut;do
  for rep in 1 2 3;do
    label=${samp}D${day}rep$rep
    mbed=$pooldir/$label.methbyread.bed.gz
    echo $mbed
    tsvs=`find $mcalldir -name "*$label*.meth.tsv"`
    subcom="cat $tsvs |\
      python $parse |\
      sort -T $pooldir -k1,1 -k2,2n |\
      bgzip > $mbed ;\
      tabix -p bed $mbed"
    echo $subcom
    eval $subcom
  done
done

