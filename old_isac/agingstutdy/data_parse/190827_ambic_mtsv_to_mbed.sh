#!/bin/bash
export LC_ALL=C
root=/mnt/d/ambic
mdir=$root/mcall
bdir=$root/mbed
[ -e $bdir ]||mkdir $bdir

parse=/home/ubuntu/Code/nanopore-methylation-utilities/mtsv2bedGraph.py

if [ "$1" == "mbed" ];then
  for mtsv in $(find $mdir -name "*tsv.gz"); do
    base=$(basename "${mtsv%.tsv.gz}")
    bases="$bases $base"
  done
  base="{}"
  in=$mdir/$base.tsv.gz
  out=$bdir/$base.bed.gz
  com="$parse -i $in | sort -S 5G -k1,1 -k2,2n -k3,3n | bgzip > $out"
  parallel $com ::: $bases
fi

if [ "$1" == "pool" ];then
  pooldir=$root/pool/mbed
  [ -e $pooldir ]||mkdir -p $pooldir
  for samp in Host StableGln StableNogln UnstableGln UnstableNogln;do
    for day in 0 30 60 90; do
      for rep in 1 2 3;do
        label=CHOZN${samp}Day${day}_$rep
        mbed=$pooldir/$label.cpg.meth.bed.gz
        echo $mbed
        beds=$(find $bdir -name "$label*bed.gz")
        gunzip -c $beds |\
          sort --parallel=94 -S 500G -k1,1 -k2,2n -k3,3n |\
          bgzip > $mbed
      done
    done
  done
fi

if [ "$1" == "bampool" ];then
  pooldir=$root/pool/bam
  [ -e $pooldir ]||mkdir -p $pooldir
  bamdir=/mnt/c/ambic/bam
  for samp in Host StableGln StableNogln UnstableGln UnstableNogln;do
    for day in 0 30 60 90; do
      for rep in 1 2 3;do
        label=CHOZN${samp}Day${day}_$rep
        bam=$pooldir/$label.pooled.bam
        echo $bam
        bams=$(find $bamdir -name "$label*bam")
        num=$(echo $bams | wc -w) 
        echo $num
        if [ $num -eq 1 ];then
          cp $bams $bam
        else
          samtools merge -@ 90 $bam $bams
        fi
        samtools index $bam
      done
    done
  done
fi

if [ "$1" == "bam_norep" ];then
  bamdir=$root/pool/bam
  outdir=/mnt/c/ambic/pooled_norep/bam
  for samp in Host StableGln StableNogln UnstableGln UnstableNogln;do
    for day in 0 30 60 90; do
      label=CHOZN${samp}Day${day}
      bam=$outdir/$label.pooled.bam
      echo $bam
      bams=$(find $bamdir -name "$label*bam")
      num=$(echo $bams | wc -w) 
      echo $num
      samtools merge -@ 90 $bam $bams
      samtools index $bam
    done
  done
fi
