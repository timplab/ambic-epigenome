#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
db=/mithril/Data/NGS/Reference/cho/ambic_sigma_IgG.fa
db=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.fa
for bam in $(find $dir -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
  fa=$dir/${base}_sigmaIgG_insertion.fasta
  out=$dir/${base}_blast.txt
  blastn -num_threads 8 -outfmt 6 \
    -db $db -query $fa -out $out
#  longrname=$(cut -f1 $out | uniq -c | sort -n | tail -n1 | awk '{ print $2 }')
#  grep $longrname $out | cut -f2,7,8,9,10
done
