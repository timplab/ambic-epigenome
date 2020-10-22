#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
for bam in $(find $dir -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
#  samtools idxstats $bam | awk '$3>0{ print $1,$3/$2*1e6 }' | sort -k2,2n | tail -n2
#  samtools view $bam | head | awk '{ print $18 }'
samtools view $bam ambic_sigma_IgG_HC ambic_sigma_IgG_LC |\
  cut -f1,2,3,4,5,8
#  awk '{ print $3,$4,$5,$(NF-1) }'
done
