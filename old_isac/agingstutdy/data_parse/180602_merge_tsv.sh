#!/bin/bash
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
pooldir=$root/pooled/methylation/meth
[ -e $pooldir ]||mkdir -p $pooldir
mcalldir=$root/mcall-cpg

for samp in Host StableGlut UnstableGlut;do
  for rep in 1 2 3;do
    label=${samp}D${day}rep$rep
    pooltsv=$pooldir/$label.meth.bed.gz
    echo $pooltsv
    tsvs=`find $mcalldir -name "*$label*.meth.tsv"`
    subcom="awk 'FNR>1{ print }' $tsvs |
    sort -T $pooldir -k1,1 -k2,2n |
    bgzip > $pooltsv; 
    tabix -p bed $pooltsv"
    echo $subcom
    eval $subcom
  done
done

