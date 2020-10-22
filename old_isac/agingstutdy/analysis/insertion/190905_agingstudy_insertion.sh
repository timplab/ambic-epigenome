#!/bin/bash
bams=$(find /mnt/*/Data -type f -name "*bam")
outdir=/mnt/d/Data/analysis
[ -e $outdir ]||mkdir $outdir

for bam in $bams; do
  base=$(basename "${bam%%.*}")
  out=$outdir/$base.supplementary.txt
  samtools view $bam ambic_sigma_IgG_HC | grep RAZU |\
    tr "\t" "\n" |\
    tr ":" "\n" |\
    tr ";" "\n" |\
    cut -d"," -f1,2,3 |\
    tr "," "\t" |\
    grep RAZU  |\
    awk '{ if ($1 != "RAZU01001824.1" || $2 < 199000 || $2 > 201000 ) print }' |\
    sort -k1,1 -k2,2n -k3,3n |\
    awk '{ print $1,int($2/100000) * 100000 }' |\
    uniq -c |\
    awk '$1>1 {print}' |\
    sort -k1,1nr > $out
done
