#!/bin/bash
export TMPDIR=/dilithium/Data/tmp
snake=181129_chosigma_snakemake.py
config=chosigmaAgingstudy_snake.config
cores=72
if [ "$1" == "makeflow" ];then
  snakemake -j $cores --snakefile $snake --configfile $config --dag |\
    dot -Tsvg > 181129_chosigma_snakemake.dag.svg
else
  snakemake -j $cores --snakefile $snake --configfile $config #-np #--unlock
fi

