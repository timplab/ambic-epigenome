#!/usr/bin/snakemake --snakefile
import pandas as pd
import os
from snakemake.utils import validate

workdir: config['workdir']
tb = pd.read_csv(config['samples'])
samples = tb.set_index("samples",drop=False)
dirs = tb.set_index("fqdir",drop=False)
fqpres = dirs.index + "/" + samples.index
samplist = list(samples.index)
dirlist = list(dirs.index)
maxthreads = 72
if cluster_config :
    print("cluster")
    maxthreads = cluster_config["__default__"]["cores"]

rule all:
    input:
        expand("bam/{sample}.sorted.bam",sample=samples.index)

rule guppy_basecall:

rule nanopolish_index:
    input:
        "{pre}.fastq.gz"
    output:
        "{pre}.fastq.gz.index.readdb"
    shell:
        "nanopolish index --verbose -d {wildcards.pre} "
        "-s {wildcards.pre}.summary.txt {input} "
        "&> {wildcards.pre}.index.log"

rule ngmlr_align:
    input:
        lambda wildcards:
            os.path.join(dirlist[samplist.index(wildcards.sample)],
            wildcards.sample+".fastq.gz")
    threads: maxthreads
    params:
        ref=config['ref']
    output:
        "bam/{sample}.sorted.bam"
    shell:
        "ngmlr -x ont --bam-fix "
        "-t {threads} -r {params.ref} -q {input} | "
        "samtools view -q 20 -b - | "
        "samtools sort -T bam/{wildcards.sample}.sorting "
        "-o {output} && "
        "samtools index {output}"

rule call_cpg:
    input:
        fq=lambda wildcards:
            os.path.join(dirlist[samplist.index(wildcards.sample)],
            wildcards.sample+".fastq.gz"),
        bam="bam/{sample}.sorted.bam"
    threads: maxthreads
    params:
        ref=config['ref']
    output:
        "mcall/{sample}.cpg.meth.tsv.gz"
    shell:
        "nanopolish call-methylation -v -t {threads} -q cpg "
        "-g {params.ref} -r {input.fq} -b {input.bam} | "
        "gzip > {output}"

rule call_gpc:
    input:
        fq=lambda wildcards:
            os.path.join(dirlist[samplist.index(wildcards.sample)],
            wildcards.sample+".fastq.gz"),
        bam="bam/{sample}.sorted.bam"
    threads: maxthreads
    params:
        ref=config['ref']
    output:
        "mcall/{sample}.gpc.meth.tsv.gz"
    shell:
        "nanopolish call-methylation -v -t {threads} -q gpc "
        "-g {params.ref} -r {input.fq} -b {input.bam} | "
        "gzip > {output}"

rule cpg_to_mbed:
    input: 
        "mcall/{sample}.cpg.meth.tsv.gz"
    output:
        "mbed/{sample}.cpg.meth.bed.gz"
    shell:
        "gunzip -c {input} | "
        "python {config[codedir]}/mtsv2bedGraph.py -m cpg --nome | "
        "sort -T tmp -k1,1 -k2,2n | "
        "bgzip > {output} && "
        "tabix -p bed {output}"

rule gpc_to_mbed:
    input: 
        "mcall/{sample}.gpc.meth.tsv.gz"
    output:
        "mbed/{sample}.gpc.meth.bed.gz"
    shell:
        "gunzip -c {input} | "
        "python {config[codedir]}/mtsv2bedGraph.py -m gpc --nome | "
        "sort -T tmp -k1,1 -k2,2n | "
        "bgzip > {output} && "
        "tabix -p bed {output}"
