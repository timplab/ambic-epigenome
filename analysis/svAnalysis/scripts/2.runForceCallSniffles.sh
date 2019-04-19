#!/bin/bash

## Sniffles Force call wrapper. Submits jobs to marcc for each .bam file in $1
## $1 = Directory containing bam files

for filename in ${1}/*.bam; do
	sbatch --partition=shared -N 1 -n 24 -t 10:0:0 -o ../slurm/`basename $filename`.slurmLog --mem=100G --mail-type=end --mail-user=kmcfarl6@jhu.edu ./forceCallSnifflesCommand.sh $filename
done
