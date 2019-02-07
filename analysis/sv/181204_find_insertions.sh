#!/bin/bash
root=/shared/Data
svdir=$root/sv
samples=$(awk 'NR>1{ print $1 }' $root/chosigma_samples.txt | tr "\n" " ")

for samp in $samples; do
  vcf=$(find $svdir -name "$samp*sniffles.vcf")
  grep IgG $vcf | cut -f1,2,5
done

