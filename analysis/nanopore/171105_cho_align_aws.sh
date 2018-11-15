#!/bin/bash
root=/shared
fqdir=${root}/fastq
bamroot=${root}/bam
align=/shared/Code/ilee/oxford/slurm/oxford_align.scr
aligner=ngmlr
#aligner=minimap2
bamdir=${bamroot}/$aligner
[ -e $bamdir ]||mkdir -p $bamdir
#ref=`readlink -f /shared/ref/picr_IgG2.fa`
ref=`readlink -f /shared/ref/picr.corrected.fa`

for fq in `find $fqdir -name "*.fq.gz"`
do
  samp=$(basename "$fq")
  samp=${samp%.fq.gz}
  echo "aligning $samp"
  prefix=${bamdir}/$samp.genomic
  log=$prefix.align.log
  echo "sbatch -o $log -e $log \
    $align -i $fq -b $prefix -a $aligner -r $ref --aws"
  sbatch -o $log -e $log \
    $align -i $fq -b $prefix -a $aligner -r $ref --aws
done

