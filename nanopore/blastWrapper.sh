#!/bin/bash
root=/atium/Data/Nanopore/oxford
wdir=/atium/Data/Nanopore/Analysis/170608_cho
base=choIgGNIH
fadir=${wdir}/fasta
dbpath1=/atium/homes/isac/Data/ambic/plasmid/4.1DHFR_VRC01WTG1M3_DGV.fa
dbpath2=/atium/homes/isac/Data/ambic/plasmid/4_0cdhfr_vrc01wtg1m3_dgv.fa
blastpath=${wdir}/blast
#mkdir $blastpath

##making database
#makeblastdb -in $dbpath2 -parse_seqids -dbtype nucl
#echo "done making database"

for runpath in `ls -d ${root}/*170530_choIgGNIH3*`
do
  samp=$(basename $runpath)
  echo $samp
  fastqpath=`ls ${runpath}/${samp}*fq.gz`
  ## convert fastq to fasta
  seqtk seq -a $fastqpath > ${fadir}/$samp.fasta

  fastapath=`ls ${fadir}/${samp}*fasta`
  echo "done making fasta"
  
  ##blasting
  echo "blasting $samp"
  blastn -num_threads 8 -outfmt 6 -db $dbpath2 -query $fastapath -out ${blastpath}/${samp}_plasmid_v2_blastn.txt
  blastn -num_threads 8 -outfmt 6 -db $dbpath1 -query $fastapath -out ${blastpath}/${samp}_plasmid_v1_blastn.txt
  echo "done blasting"
done
