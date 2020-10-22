#!/bin/bash
root=/kyber/Data/NGS/projects/ambic/agingstudy/atacseq/190211_choSigmaD0_atac/data/atacseq
gs=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.fa.fai
beddir=$root/bed
for bed in $(find $beddir -name "*bed.gz"); do
  base=$(basename "$bed")
  base=${base%%.*}
  bases="$bases $base"
  out=$root/genomecov/$base.genomecov.bed
done
  out=$root/genomecov/{}.genomecov.bed
  com="bedtools genomecov -i $beddir/{}.bed.gz -g $gs > $out"
  parallel $com ::: $bases
