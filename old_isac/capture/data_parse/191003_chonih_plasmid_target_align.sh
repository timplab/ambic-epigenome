#!/bin/bash
dir=/home/isac/Data/ambic/targeted/191001_CHO_NIH_IgG3_targeted
base=191001_CHO_NIH_IgG3_targeted

if [ "$1" == "align" ];then
  idx=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.nihIgG.mmi
  for fq in $(find $dir -name "*fastq"); do
    base=$(basename "$fq")
    base=${base%%.*}
    echo $base
    out=$dir/$base.sorted.bam
    log=$dir/$base.align.log
    minimap2 --MD -t 64 -ax map-ont -L $idx $fq 2> $log |\
      samtools view -b -q 20 - |\
      samtools sort -o $out
    samtools index $out
    samtools idxstats $out
  done
fi

if [ "$1" == "getontarg" ];then
  bam=$dir/$base.sorted.bam
  plas=4_0cdhfr_vrc01wtg1m3_dgv 
  samtools view $bam $plas | grep RAZU | awk '{ print $(NF-2) }' |\
    cut -d":" -f3 | tr ";" "\n" | grep . | awk 'BEGIN{ FS=",";OFS="\t"}{ print $1,$2; }' |\
    grep -v $plas  > $dir/insert_loc.txt
fi
