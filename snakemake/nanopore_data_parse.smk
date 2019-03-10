#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate

rule ngmlr_align:
	input:
		"{dir}/reads/{sample}.fastq.gz"
	params:
		config['reference']
	threads: maxthreads
	output:
		"{dir}/bam/{sample}.sorted.bam"
	shell:
		"ngmlr -x ont --bam-fix "
		"-t {threads} -r {params} -q {input} | "
		"samtools view -q 20 -b - | "
		"samtools sort -T {wildcards.dir}/bam/{wildcards.sample}.sorting "
		"-o {output} && "
		"samtools index {output}"

rule nanopolish_index:
	input: 
		"{pre}.fastq.gz" 
	output: 
		"{pre}.fastq.gz.index.readdb" 
	shell:
		"nanopolish index --verbose -d {wildcards.pre} "
#		"-s {wildcards.pre}.summary.txt "  # use this if summary provided
		"{input} "
		"&> {wildcards.pre}.index.log"

rule call_cpg:
	input:
		fq="{dir}/reads/{sample}.fastq.gz",
		bam="{dir}/bam/{sample}.sorted.bam",
		db="{dir}/reads/{sample}.fastq.gz.index.readdb"
	params:
		config['reference']
	threads: maxthreads
	output:
		"{dir}/mcall/{sample}.cpg.meth.tsv.gz"
	shell:
		"nanopolish call-methylation -v -t {threads} -q cpg "
		"-g {params} -r {input.fq} -b {input.bam} | " 
		"gzip > {output}"

rule mtsv_to_mbed:
	input:
		"{dir}/mcall/{sample}.cpg.meth.tsv.gz"
	params:
		config['codedir']
	output:
		"{dir}/mbed/{sample}.cpg.meth.bed.gz"
	shell:
		"gunzip -c {input} | "
		"python {params}/scripts/mtsv2bedGraph.py | "
		"sort -T tmp -k1,1 -k2,2n | "
		"bgzip > {output} && "
		"tabix -p bed {output}"

rule mbed_to_mfreq: 
	input: 
		"{dir}/mbed/{sample}.{mod}.meth.bed.gz" 
	params: 
		config['codedir'] 
	output: 
		mfreq="{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz", 
		log="{dir}/mfreq/{sample}.{mod}.mfreq.log" 
	shell: 
		"python -u {params}/scripts/parseMethylbed.py frequency -v " 
		"-i {input} -m {wildcards.mod} 2> {output.log} | " 
		"bgzip > {output.mfreq} && " 
		"tabix -b 2 -e 2 {output.mfreq}"

rule mfreq_to_wig: 
	input: 
		"{dir}/mfreq/{sample}.{mod}.mfreq.txt.gz" 
	params: 
		config['codedir'] 
	output: 
		methwig=temp("{dir}/bigwig/{sample}.{mod}.methylation.wig"), 
		covwig=temp("{dir}/bigwig/{sample}.{mod}.methcoverage.wig"), 
		log="{dir}/bigwig/{sample}.{mod}.wig.log" 
	shell: 
		"python {params}/script/makeWig.py -v -i {input} " 
		"-o {output.methwig} -c {output.covwig} &> {output.log}"

rule wig_to_bigwig: 
	input: 
		"{dir}/bigwig/{sample}.{mod}.{type}.wig", 
		"data/hg38/hg38_genomesize.txt" 
	output: 
		"{dir}/bigwig/{sample}.{mod}.{type}.bw" 
	shell: 
		"wigToBigWig {input} {output}"
