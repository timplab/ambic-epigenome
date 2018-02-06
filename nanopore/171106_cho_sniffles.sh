#!/bin/bash
root=/shared
bamdir=${root}/mergedbam
svdir=${root}/sniffles
ref=/shared/ref/picr_IgG2.fa
slurmscr=/shared/Code/ilee/sv/slurm/sniffles.scr
srcdir=/shared/Code

for bam in `find $bamdir -name "*.pooled.bam"`
do
  samp=$(basename "$bam")
  samp=${samp%%.*}
  log=$svdir/$samp.sniffles.log
  echo "sbatch --oversubscribe $slurmscr -s $srcdir -b $bam -p $samp -o $svdir --xargs "-s 3 -n 10""
  sbatch --oversubscribe --error=$log --output=$log $slurmscr -s $srcdir -b $bam -p $samp -o $svdir --xargs "-s 3 -n 10"
done

