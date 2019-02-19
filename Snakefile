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
include:
	smkdir + "atacseq_data_parse.smk"

nanopore_tb = pd.read_csv(config['codedir']+"/nanopore_sample_info.csv")
atacseq_tb = pd.read_csv(config['codedir']+"/atacseq_sample_info.csv")

rule parse_nanopore:
  input:
    expand("data/nanopore/bam/{sample}.sorted.bam",
			sample=nanopore_tb['sample'])

rule parse_atacseq:
	input:
		"data/atacseq/peaks/individual_peaks_counts.txt"
#		expand("data/atacseq/bam/{sample}.nodup.sorted.bam",
#			sample=atacseq_tb['sample'])
		
