#!/bin/bash
assembly=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_singler_2line.fa
dir=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly
ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
out=$dir/200205_homology_blastn.txt

if [ ! -e $ref.nsd ];then
  makeblastdb -in $ref -dbtype nucl -parse_seqids
fi

blastn -num_threads 8 -outfmt 6 \
  -db $ref -query $assembly -out $out
