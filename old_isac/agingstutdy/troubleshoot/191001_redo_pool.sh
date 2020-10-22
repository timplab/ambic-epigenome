#!/bin/bash
export LC_ALL=C
root=/mnt/d/Data/redo
outroot=$root/pooled_norep
indir=$root/pooled_rep

if [ "$1" == "mbed" ];then
  outdir=$outroot/mbed
  [ -e $outdir ]||mkdir $outdir
  bases=$(find $indir/mbed -name "*bed.gz" | tr "/" " " | awk '{print $NF}' | tr "_" "\t" | cut -f1 | sort | uniq)
  base="{}"
  inbed=$indir/mbed/$base*bed.gz
  out=$outdir/$base.cpg.meth.bed.gz
  tmpdir=$root/tmp/mbed/$base
  com="mkdir -p $tmpdir &&\
    gunzip -c $inbed | sort -T $tmpdir -k1,1 -k2,2n -k3,3n | bgzip > $out &&\
    tabix -p bed $out"
#  com="echo $base"
  parallel $com ::: $bases
fi

if [ "$1" == "mfreq" ];then
  scr="/home/ubuntu/Code/nanopore-methylation-utilities/parseMethylbed.py frequency"
  outdir=$outroot/mfreq
  [ -e $outdir ]||mkdir $outdir
  bases=$(find $outroot/mbed -name "*bed.gz" | tr "/" " " | awk '{print $NF}' | tr "." "\t" | cut -f1)
  base="{}"
  inbed=$outroot/mbed/$base.cpg.meth.bed.gz
  out=$outdir/$base.cpg.meth.mfreq.txt.gz
  com="gunzip -c $inbed | $scr | bgzip > $out"
#  com="echo $base"
  parallel $com ::: $bases
fi

if [ "$1" == "bam" ];then
  outdir=$outroot/bam
  [ -e $outdir ]||mkdir $outdir
  bases=$(find $indir/bam -type f -name "*bam" | tr "/" " " | awk '{print $NF}' | tr "_" "\t" | cut -f1 | sort | uniq)
  for base in $bases; do 
    out=$outdir/$base.pooled.bam
    inbam=$indir/bam/$base*bam
    com="samtools merge -@ 96 $out $inbam &&\
      samtools index -@ 96 $out"
    echo $com
    eval $com
  done
fi

if [ "$1" == "bed" ];then
  # first rep
  outdir=$indir/bed
  [ -e $outdir ]||mkdir $outdir
  bases=$(find $root/../bed -type f -name "*bed" | tr "/" " " | awk '{print $NF}' | tr "_" "\t" | awk '{print $1"_"$2}' | sort | uniq)
  base="{}"
  tmpdir=$root/tmp/bed/$base
  inbed=$root/../bed/$base*bed
  out=$outdir/$base.pooled.bed.gz
  com="mkdir -p $tmpdir &&\
    cat $inbed | sort -T $tmpdir -k1,1 -k2,2n -k3,3n | bgzip > $out"
#  com="echo $base"
  echo "rep"
  parallel $com ::: $bases
  # second norep
  outdir=$outroot/bed
  [ -e $outdir ]||mkdir $outdir
  bases=$(find $indir/bed -type f -name "*bed.gz" | tr "/" " " | awk '{print $NF}' | tr "_" "\t" | cut -f1 | sort | uniq)
  base="{}"
  tmpdir=$root/tmp/bed/$base
  inbed=$indir/bed/$base*bed.gz
  out=$outdir/$base.pooled.bed.gz
  com="mkdir -p $tmpdir &&\
    gunzip -c $inbed | sort -T $tmpdir -k1,1 -k2,2n -k3,3n | bgzip > $out"
#  com="echo $base"
  echo "norep"
  parallel $com ::: $bases
fi

if [ "$1" == "bincov" ];then
  outdir=$outroot/binned_coverage
  [ -e $outdir ]||mkdir $outdir
  bins=$outdir/picr_sigma_bins.bed
  if [ ! -e $bins ];then
    binner=/home/ubuntu/Code/ambic-epigenome/scripts/bin_genome.py
    ref=/mnt/d/Data/Reference/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.sigmaIgG.fa
    $binner -f $ref.fai -b 10000 > $bins
  fi
  bases=$(find $outroot/bed -type f -name "*bed.gz" | tr "/" " " | awk '{print $NF}' | tr "." "\t" | cut -f1 | sort | uniq)
  base="{}"
  bed=$outroot/bed/$base.pooled.bed.gz
  cov=$outdir/$base.binned_coverage.bed
  com="bedtools coverage -a $bins -b $bed > $cov"
  echo $com
  parallel $com ::: $bases
fi



