#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate 

workdir: config['workdir']
maxthreads = 72 
samples = pd.read_table(config['samples']).set_index("samples",drop=False)
if config['suffix'] :
    samples.index = [ x+"."+config['suffix'] for x in samples.index ]

rule all:
    input: 
#        "mfreq/"+config['base']+".meth.bsmooth.rds"
        expand("mfreq_pooled/{sample}.methfreq.txt.gz",sample=samples.index)

rule pool_to_replicates:
    params:
        codedir=config['codedir']
    output:
        "mcall_pooled/{sample}.meth.tsv"
    shell:
        "{params.codedir}/util/merge_methylation_tsv.sh "
        "-d mcall_replicates -b {wildcards.sample} > {output}"

rule convert_mtsv_to_mbed:
    input:
        "mcall_pooled/{sample}.meth.tsv"
    params:
        codedir = config['codedir']
    output:
        "mbed_pooled/{sample}.meth.bed.gz"
    shell:
        "python {params.codedir}/methylation/mtsv2bedGraph_nostrand.py "
        "-i {input} | sort -k1,1 -k2,2n | bgzip > "
        "{output} && tabix -p bed {output}"

rule get_methylation_frequency:
    input:
        "mbed_pooled/{sample}.meth.bed.gz"
    params:
        codedir = config['codedir']
    output:
        "mfreq_pooled/{sample}.methfreq.txt.gz"
    log:
        "log/mfreq/{sample}.mfreq.log"
    shell:
        "python {params.codedir}/methylation/parseMethylbed.py "
        "frequency --verbose -i {input} 2> {log} | "
        "bgzip > {output} && tabix -b 2 -e 2 {output}"

#rule smooth_methylation:
#    input:
#        expand("mfreq_pooled/{sample}.methfreq.txt",sample=samples.index)
#    params:
#        codedir = config['codedir'],
#        pdpath = config['samples']
#    output:
#        "mfreq/"+config['base']+".meth.bsmooth.rds"
#    shell:
#        "{params.codedir}/methylation/bsmooth_mfreq.R "
#        "-p {params.pdpath} -o {output} {input}"
