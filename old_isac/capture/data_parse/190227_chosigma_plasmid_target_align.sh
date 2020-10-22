#!/bin/bash
dir=/home/isac/Data/ambic/targeted/190225_chosigma_targeted
base=CHOZNStableNoglnDay0rep3
fq=$dir/$base.fastq.gz
bam=$dir/$base.sorted.bam

if [ "$1" == "align" ];then
  idx=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.sigmaIgG.mmi
  log=$dir/$base.align.log
  minimap2 --MD -t 10 -ax map-ont -L $idx $fq 2> $log |\
    samtools view -b -q 20 - |\
    samtools sort -o $bam
  samtools index $bam
fi

if [ "$1" == "getontarg" ];then
  bam=$dir/$base.sorted.bam
  plas=ambic_sigma_IgG_HC
  samtools view $bam $plas | grep RAZU | awk '{ print $(NF-2) }' |\
    cut -d":" -f3 | tr ";" "\n" | grep . | awk 'BEGIN{ FS=",";OFS="\t"}{ print $1,$2; }' |\
    grep -v sigma  > $dir/insert_loc.txt
fi
