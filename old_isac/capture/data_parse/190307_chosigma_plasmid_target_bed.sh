#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
cd $dir
for bam in $(find * -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
  bed=$dir/$base.bed
  bedtools bamtobed -i $bam > $bed
done
