#!/bin/bash
root=/shared
codedir=$root/Code
repo=$codedir/ilee
srcpath=$repo/dnamods
slurmpath=$repo/oxford/slurm
logroot=$root/log
rawtardir=$root/rawtar
rawdir=$root/raw
unorgdir=$root/raw_unorg

for tag in host IgG; do
  for rep in 1 2 3; do
    wdir=$unorgdir/$tag/$rep
    n=10000
    base=2017
    out=$rawdir/$tag/$rep
    [ -e $out ]||mkdir -p $out
    $repo/oxford/combineRuns.py -n $n -s $base -o $out $wdir
  done
done
