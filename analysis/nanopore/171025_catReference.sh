#!/bin/bash

ref=/mithril/Data/NGS/Reference/cho/criGri1.fa
plasmid=/atium/homes/isac/Data/ambic/plasmid/4_0cdhfr_vrc01wtg1m3_dgv.fa
version=2

dir=$(dirname "$ref")
base=$(basename "$ref")
base=${base%.fa}

out=${dir}/${base}_plasmid$version.fa
cat $ref $plasmid > $out

echo "finished concatenating into $out"


