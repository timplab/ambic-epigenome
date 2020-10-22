#!/bin/bash
bams=$(find /mnt/*/Data -type f -name "*bam")

for bam in $bams; do
  root=${bam%/pooled_norep*}
  base=$(basename "${bam%%.*}")
  outdir=$root/sniffles
  [ -e $outdir ]||mkdir $outdir
  vcf=$outdir/$base.sniffles.vcf
  if [ -e $vcf ];then
    continue
  fi
  log=$outdir/$base.sniffles.log
  sniffles -m $bam -v $vcf --tmp_file $outdir \
    -s 5 -t 48 -l 50 -n -1 --genotype --cluster &> $log
done
