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

rule nanopore_sv:
  input:
	   vcf="data/nanopore/sv/merged_final_SURVIVOR_1kbpdist_typesave.vcf"
#    expand("data/nanopore/sv/{sample}.minimap2.SURVIVOR.sorted.vcf",
#			sample=nanopore_tb['sample'])

rule nanopore_methylation:
  input:
    expand("data/nanopore/mfreq/{sample}.cpg.mfreq.txt.gz",
			sample=nanopore_tb['sample'])

rule parse_atacseq:
	input:
		out="data/atacseq/peaks/individual_peaks_counts.txt"
#		expand("data/atacseq/bam/{sample}.nodup.sorted.bam",
#			sample=atacseq_tb['sample'])
		
