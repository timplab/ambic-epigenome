#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
cd $dir
for bam in $(find * -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
  vcf=$dir/$base.sniffles.vcf
  log=$dir/$base.sniffles.log
  sniffles -m $bam -v $vcf \
    --tmp_file $dir/$base.tmp \
    -s 10 -n -1 --genotype --cluster &> $log
done
