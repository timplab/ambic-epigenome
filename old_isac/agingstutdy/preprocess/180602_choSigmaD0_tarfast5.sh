#!/bin/bash
root=/shared/data/raw/fast5
day=0

batchscr="/shared/Code/ilee/slurm/batchcommand.scr"

for samp in Host StableGlut UnstableGlut;do
  for rep in 1 2 3;do
    label=${samp}D${day}rep$rep
    echo $label
    mkdir $root/$label
    labdir=`find $root -maxdepth 1 -name "*$label*"`
    mv $labdir $root/$label
    cd $root
    tar -vczf $root/$label.fast5.tgz $label > $root/$label.fast5list.txt
  done
done
