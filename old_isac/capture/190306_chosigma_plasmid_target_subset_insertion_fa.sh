#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
for bam in $(find $dir -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  echo $base
  vcf=$dir/$base.sniffles.vcf
  qnames=$dir/$base.insert.qnames.txt
  samtools view -h $bam ambic_sigma_IgG_HC ambic_sigma_IgG_LC | cut -f1 > $qnames
  grep sigma $vcf | cut -d";" -f8 |\
    sed "s/RNAMES=//" |\
    tr "," "\n" | sort >> $qnames
  insertsam=$dir/$base.insert.sam
  sam=$dir/$base.sam
  fa=$dir/${base}_sigmaIgG_insertion.fasta
  samtools view -H $bam > $insertsam
  samtools view $bam > $sam 
  awk 'NR==FNR{c[$1]++;next};c[$1]' $qnames $sam >> $insertsam
  samtools fasta $insertsam > $fa
done
