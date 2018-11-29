#!/bin/bash
snake=181129_chosigma_snakemake.py
config=chosigmaAgingstudy_snake.config
if [ "$1" == "makeflow" ];then
  snakemake -j 8 --snakefile $snake --configfile $config --dag |\
    dot -Tsvg > 181129_chosigma_snakemake.dag.svg
else
  snakemake -j 8 --snakefile $snake --configfile $config  #-np #--unlock
fi

