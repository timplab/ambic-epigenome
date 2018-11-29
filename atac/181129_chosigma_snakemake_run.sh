#!/bin/bash
snake=181129_chosigma_snakemake.py
config=chosigmaAgingstudy_snake.config
snakemake -j 8 --snakefile $snake --configfile $config 
