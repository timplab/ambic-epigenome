#!/bin/bash
root=/scratch/users/ilee29@jhu.edu/171025_cho
fqdir=${root}/fastq
bamdir=${root}/bam
[ -d $bamdir ] || mkdir $bamdir

align=/home-2/ilee29@jhu.edu/Code/ilee/oxford/slurm/oxford_align.scr
aligner=ngmlr
ref=/scratch/groups/wtimp1/Reference/cho/criGri1_plasmid2.fa

for fq in `find $fqdir -name "*fq.gz"`
do
  samp=$(basename "$fq")
  samp=${samp%.fq.gz}
  echo "aligning $samp"
  log=$samp.align.log
  sbatch -D $bamdir --time=24:0:0 --output $log --error $log \
    $align -i $fq -a $aligner -r $ref --marcc
done

