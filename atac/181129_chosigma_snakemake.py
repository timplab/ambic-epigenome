#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate 

workdir: config['workdir']
samples = pd.read_table(config['samples']).set_index("samples",drop=False)

rule all:
    input: 
        expand("fastq_trimmed/{sample}_trimmed.fq.gz",sample=samples.index)

rule trim:
    input :
        "fastq/{sample}.fastq.gz"
    output :
        "fastq_trimmed/{sample}_trimmed.fq.gz"
    log:
        "log/trim/{sample}.trim.log"
    shell :
        "trim_galore {input} -o fastq_trimmed &> {log}"
