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
        "coverage/"+config['base']+".coverage.bed"
#        "sv/"+config['base']+".forcecall.merged.vcf"
#@        expand("sv/{sample}.sniffles.forcecall.vcf",sample=samples.index)

# sniffles
rule sniffles_detect_svs:
    input:
        "bam/{sample}.bam"
    output:
        "sv/{sample}.sniffles.vcf"
    log:
        "sv/{sample}.sniffles.log"
    shell:
        "sniffles -m {input} -v {output} "
        "--tmp_file sv/{wildcards.sample}.tmp "
        "-s 2 -n -1 --genotype --cluster &> {log}"

rule merge_sniffles_across_samples:
    input:
        expand("sv/{sample}.sniffles.vcf",sample=samples.index)
    output:
        fofn = "sv/"+config['base']+".sniffles.fofn",
        vcf = "sv/"+config['base']+".merged.vcf"
    log:
        "log/sv/"+config['base']+".sv.firstmerge.log"
    shell:
        "echo {input} | tr ' ' '\n' > {output.fofn} && "
        "SURVIVOR merge {output.fofn} 1000 1 1 -1 -1 -1 "
        "{output.vcf} &> {log}"

rule sniffles_force_detect_svs:
    input:
        vcf = "sv/"+config['base']+".merged.vcf",
        bam = "bam/{sample}.bam"
    output:
        "sv/{sample}.sniffles.forcecall.vcf"
    log:
        "log/sv/{sample}.forcecall.log"
    shell:
        "sniffles -m {input.bam} -v {output} "
        "--tmp_file sv/{wildcards.sample}.tmp "
        "--Ivcf {input.vcf} -s 2 -n -1 --genotype "
        "--cluster &> {log}"
    
rule merge_forecall_sniffles_across_samples:
    input:
        expand("sv/{sample}.sniffles.forcecall.vcf",sample=samples.index)
    output:
        fofn = "sv/"+config['base']+".sniffles.forcecall.fofn",
        vcf = "sv/"+config['base']+".forcecall.merged.vcf"
    log:
        "log/sv/"+config['base']+".sv.secondmerge.log"
    shell:
        "echo {input} | tr ' ' '\n' > {output.fofn} && "
        "SURVIVOR merge {output.fofn} 1000 -1 1 -1 -1 -1 "
        "{output.vcf} &> {log}"


# coverage analysis
rule bam_to_bed:
    input:
        "bam/{sample}.bam"
    output:
        "bed/{sample}.bed"
    shell:
        "bedtools bamtobed -i {input} | "
        "sort -k1,1 -k2,2n > {output}"

rule bin_genome:
    input:
        config['genome']
    params:
        codedir = config['codedir']
    output:
        config['genome']+".bins.bed"
    shell:
        "{params.codedir}/script/bin_genome.sh {input} {output} 10000"

rule calculate_coverage:
    input:
        reads = "bed/{sample}.bed",
        bins = config['genome']+".bins.bed"
    output:
        "coverage/{sample}.coverage.bed"
    shell:
        "bedtools coverage -a {input.bins} -b {input.reads} > {output}"

rule merge_coverage:
    input:
        expand("coverage/{sample}.coverage.bed",sample=samples.index)
    params:
        codedir = config['codedir']
    output:
        "coverage/"+config['base']+".coverage.bed"
    shell:
        "{params.codedir}/util/make_coverage_matrix_from_bedtools_coverage.sh {input} > {output}"
