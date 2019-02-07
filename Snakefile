#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

configfile:
  "snakemake_config.yml"
workdir:
  config['workdir']
maxthreads = config['threads']
smkdir = config['codedir'] + "/snakemake/"
include:
  smkdir + "nanopore_data_parse.smk"

nanopore_tb = pd.read_csv(config['codedir']+"/nanopore_sample_info.csv")

rule parse_nanopore:
  input:
    expand("data/nanopore/bam/{sample}.sorted.bam",
			sample=nanopore_tb['sample'])
