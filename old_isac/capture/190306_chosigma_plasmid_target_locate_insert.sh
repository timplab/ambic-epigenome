#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
for bam in $(find $dir -name "CHOZNUnstableNoglnDay0rep3*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
  python insert_locater.py -b $bam -c ambic_sigma_IgG_HC ambic_sigma_IgG_LC
done
