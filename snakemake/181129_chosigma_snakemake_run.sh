#!/bin/bash
tmpdir=/shared/Data/tmp
[ -e $tmpdir ]||mkdir -p $tmpdir
export TMPDIR=/shared/Data/tmp
export LC_ALL=C
cores=72
# arg options : atacseq, genome
snake=181129_chosigma_${1}_snakemake.py
config=chosigmaAgingstudy_${1}_snakemake.config
dag=181129_chosigma_${1}.dag.svg
snakemake -j $cores --snakefile $snake --configfile $config #-np #--unlock
#snakemake -j $cores --snakefile $snake --configfile $config --dag |\
#  dot -Tsvg > $dag
