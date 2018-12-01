#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate 

workdir: config['workdir']
maxthreads = 72 
samples = pd.read_table(config['samples']).set_index("samples",drop=False)

rule all:
    input: 
        "peaks/"+config['base']+"_peaks_counts.txt"
        
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

rule filter_reads_with_multiple_alignments:
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
                samtools index {output} && \
                rm bam/{wildcards.sample}.dupmark.bam"
rule bam_to_bed:
    input:
        "bam/{sample}.nodup.sorted.bam"
    params:
        codedir=config['codedir']
    output:
        "bed/{sample}.bed.gz"
    shell:
        "{params.codedir}/script/atacseq_bam_to_bed.sh -i {input} -o {output}"

rule merge_bed:
    input:
        expand("bed/{sample}.bed.gz",sample=samples.index)
    threads: maxthreads
    output:
        "bed/"+config['base']+".merged.bed.gz"
    shell:
        "export LC_ALL=C;\
                gunzip -c {threads} {input} |\
                sort -k1,1 -k2,2n |\
                gzip > {output}"

rule find_peaks_from_pool:
    input:
        "bed/"+config['base']+".merged.bed.gz"
    params:
        base=config['base'],
        codedir=config['codedir']
    output:
        "peaks/"+config['base']+"_peaks.saf"
    log:
        "log/"+config['base']+".peaks.log"
    shell:
        "macs2 callpeak -t {input} -f BED -n peaks/{params.base} "
        "-g mm -p 0.01 --shift -75 --extsize 150 "
        "--nomodel -B --SPMR --keep-dup all --call-summits &> {log} && "
        "{params.codedir}/script/atacseq_bed_to_saf.sh "
        "peaks/{params.base}_peaks.narrowPeak > {output}"

rule get_feature_counts:
    input:
        reads=expand("bam/{sample}.nodup.sorted.bam",sample=samples.index),
        peaks="peaks/"+config['base']+"_peaks.saf"
    threads: 64
    output:
        "peaks/"+config['base']+"_peaks_counts.txt"
    log:
        "log/"+config['base']+".featurecounts.log"
    shell:
        "featureCounts --verbose -a {input.peaks}"
        "-o {output} -F SAF -T {threads} {input.reads} &> {log}"
