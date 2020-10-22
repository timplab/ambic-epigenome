#!/bin/bash
root=/mnt/e
fqdir=$root/fastq
bamdir=$root/bam
rdir=$root/reads
mdir=$root/mcall
mi=$root/Reference/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.sigmaIgG.mmi
ref=$root/Reference/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.sigmaIgG.fa
[ -e $bamdir ]||mkdir $bamdir

if [ "$1" == "align" ];then
  for fq in $(find $fqdir -name "*fastq.gz" ! -name "*Day60*" ! -name "*Day90*"); do
    base=$(basename "${fq%%.*}")
    bam=$bamdir/$base.sorted.bam
    log=$bamdir/$base.align.log
    if [ ! -e $bam ];then
      echo $base
      minimap2 --MD -I 16G -t 90 -a $mi $fq 2> $log |\
        samtools view -bT $ref - |\
        samtools sort -@ 6 -o $bam &&
        samtools index $bam
    else 
      echo $base already aligned
    fi
  done
fi

for dir in $(find $rdir/* -maxdepth 0 -type d); do
  samps="$samps $(basename "$dir")"
done

np=~/Code/nanopolish/nanopolish
[ -e $mdir ]||mkdir $mdir
for samp in $samps; do
  fq=$fqdir/$samp.fastq.gz
  db=$fq.index.readdb
  if [ ! -e $db ];then
    index="$index $samp"
  else
    echo $samp already indexed
  fi
done
samp="{}"
fq=$fqdir/$samp.fastq.gz
dir=$rdir/$samp
com="$np index -d $dir $fq"
parallel $com ::: $index

for samp in $samps; do
  fq=$fqdir/$samp.fastq.gz
  db=$fq.index.readdb
  echo $samp
  bam=$bamdir/$samp.sorted.bam
  tsv=$mdir/$samp.cpg.meth.tsv.gz
  log=$mdir/$samp.cpg.mcall.log
  $np call-methylation -v -t 96 \
    -r $fq -b $bam -g $ref 2> $log |\
    gzip > $tsv
  
done
