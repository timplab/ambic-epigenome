#!/bin/bash
# blast the sigma plasmid to the chok1 genome
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis
insertpath=$root/insertion
[ -e $insertpath ]||mkdir $insertpath

faroot=/mithril/Data/NGS/Reference/cho
if [ "$1" == "plasmid2genome" ];then
  query=$faroot/ambic_sigma_IgG.fa
  db=$faroot/chok1/genome/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa
  blast=$insertpath/ambicSimgaIgG_to_chok1_blastn.txt
elif [ "$1" == "genome2plasmid" ];then
  db=$faroot/ambic_sigma_IgG.fa
  query=$faroot/chok1/genome/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa
  blast=$insertpath/chok1_to_ambicSimgaIgGblastn.txt
  arg="-num_threads 8"
fi


if [ ! -e $db.nsq ];then
  echo "making blast database"
  makeblastdb -in $db -parse_seqids -dbtype nucl
fi

echo "blasting"
blastn $arg -outfmt 6 -db $db -query $query -out $blast
