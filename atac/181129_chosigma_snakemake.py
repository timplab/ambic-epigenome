#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate 

workdir: config['workdir']
threads: config['threads']
maxthreads = 36
samples = pd.read_table(config['samples']).set_index("samples",drop=False)

rule all:
    input: 
        expand("bam/{sample}.nodup.sorted.bam",sample=samples.index)
        
rule trim:
    input:
        "fastq/{sample}.fastq.gz"
    output:
        "fastq_trimmed/{sample}_trimmed.fq.gz"
    log:
        "log/trim/{sample}.trim.log"
    shell:
        "trim_galore {input} -o fastq_trimmed &> {log}"

rule align:
    input:
        "fastq_trimmed/{sample}_trimmed.fq.gz"
    params:
        btidx=config['bowtieidx']
    threads: maxthreads # max 36, controlled by "-j" option into command or "threads" config parameter
    output:
        temp("bam/{sample}.bam")
    log:
        "log/align/{sample}.align.log"
    shell:
        "bowtie2 -k 4 -p {threads} -t --local \
                -x {params.btidx} -U {input} 2> {log} |\
                samtools view -Sb - > {output}"

rule filter_multimap:
    input:
        "bam/{sample}.bam"
    params:
        codedir=config['codedir']
    output:
        temp("bam/{sample}.multimapfiltered.bam")
    shell:
        "samtools sort -n -T bam/{wildcards.sample}.qsorting {input} |\
                samtools view -h - |\
                python {params.codedir}/script/assign_multimappers.py -k 4 |\
                samtools view -F 1804 -Su - |\
                samtools sort -T bam/{wildcards.sample}.mpsorting -o {output}"

rule remove_duplicates:
    input:
        "bam/{sample}.multimapfiltered.bam"
    threads: maxthreads
    output:
        "bam/{sample}.nodup.sorted.bam"
    log:
        "log/removedup/{sample}.dupremoval.log"
    shell:
        "picard -Xmx8G -Xms256M \
                -XX:ParallelGCThreads={threads} \
                -Djava.io.tmpdir=bam/dupremove \
                MarkDuplicates \
                INPUT={input} OUTPUT=bam/{wildcards.sample}.dupmark.bam \
                METRICS_FILE=bam/{wildcards.sample}.duplicate_metrics.txt \
                VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true REMOVE_DUPLICATES=false \
                &> {log} &&\
                samtools view -F 1804 -b bam/{wildcards.sample}.dupmark.bam > {output} &&\
                samtools index {output}"

