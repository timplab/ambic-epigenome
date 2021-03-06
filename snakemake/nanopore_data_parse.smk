#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate
nanopore_tb = pd.read_csv(config['codedir']+"/nanopore_sample_info.csv",comment="#")

##################################################
#
# alignment
#
##################################################

rule ngmlr_align:
	input:
		"reads/{sample}.fastq.gz"
	params:
		config['reference']
	threads: maxthreads
	output:
		"bam/{sample}.ngmlr.sorted.bam"
	log:
		"bam/{sample}.ngmlr.align.log"
	shell:
		"ngmlr -x ont --bam-fix "
		"-t {threads} -r {params} -q {input} 2> {log} | "
		"samtools view -q 20 -b - | "
		"samtools sort -T bam/{wildcards.sample}.sorting "
		"-o {output} && "
		"samtools index {output}"
	
rule minimap2_makeindex:
	input:
		config['reference']
	output:
		os.path.splitext(config['reference'])[0] + ".map_ont.mmi"
	shell:
		"minimap2 -x map-ont -d {output} {input}"

rule minimap2_align:
	input:
		os.path.splitext(config['reference'])[0] + ".map_ont.mmi",
		"reads/{sample}.fastq.gz"
	threads: maxthreads
	output:
		"bam/{sample}.minimap2.sorted.bam"
	log:
		"bam/{sample}.minimap2.align.log"
	shell:
		"minimap2 --MD -L -t {threads} -a {input} 2> {log} | "
		"samtools view -q 20 -b - | "
		"samtools sort -T bam/{wildcards.sample}.sorting "
		"-o {output} && "
		"samtools index {output}"

##################################################
#
# methylation calling
#
##################################################

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
		fq="reads/{sample}.fastq.gz",
		bam="bam/{sample}.minimap2.sorted.bam", # change "minimap2"/"ngmlr" to change aligner
		db="reads/{sample}.fastq.gz.index.readdb"
	params:
		config['reference']
	threads: maxthreads
	output:
		"methylation/{sample}.cpg.meth.tsv.gz"
	shell:
		"nanopolish call-methylation -v -t {threads} -q cpg "
		"-g {params} -r {input.fq} -b {input.bam} | " 
		"gzip > {output}"

rule mtsv_to_mbed:
	input:
		meth="methylation/{sample}.cpg.meth.tsv.gz",
		fa=config['reference']
	params:
		config['codedir']
	output:
		"methylation/{sample}.cpg.meth.bed.gz"
	shell:
		"python {params}/nanopore-methylation-utilities/mtsv2bedGraph.py -g {input.fa} -i {input.meth} | "
		"sort -T tmp -k1,1 -k2,2n | "
		"bgzip > {output} && "
		"tabix -p bed {output}"

rule bam_methylation:
	input:
		meth="methylation/{sample}.cpg.meth.bed.gz",
		bam="bam/{sample}.minimap2.sorted.bam"
	threads:
		maxthreads
	params:
		config['codedir']
	output:
		"methylation/{sample}.cpg.meth.bam"
	shell:
		"python {params}/nanopore-methylation-utilities/convert_bam_for_methylation.py "
        "-t {threads} -c {input.meth} -b {input.bam} | "
		"samtools sort -o {output} && "
		"samtools index {output}"
        
rule mbed_to_mfreq: 
	input: 
		"methylation/{sample}.{mod}.meth.bed.gz" 
	params: 
		config['codedir'] 
	log: 
		"methylation/{sample}.{mod}.mfreq.log" 
	output: 
		"methylation/{sample}.{mod}.mfreq.txt.gz"
	shell: 
		"python -u {params}/nanopore-methylation-utilities/parseMethylbed.py frequency -v " 
		"-i {input} -m {wildcards.mod} 2> {log} | " 
		"bgzip > {output} && " 
		"tabix -b 2 -e 2 {output}"

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
		"python {params}/scripts/makeWig.py -v -i {input} " 
		"-o {output.methwig} -c {output.covwig} &> {output.log}"

rule wig_to_bigwig: 
	input: 
		"{dir}/bigwig/{sample}.{mod}.{type}.wig", 
		"data/hg38/hg38_genomesize.txt" 
	output: 
		"{dir}/bigwig/{sample}.{mod}.{type}.bw" 
	shell: 
		"wigToBigWig {input} {output}"

##################################################
#
# SV detection 
#
##################################################

rule sniffles_sv:
	input:
		"bam/{sample}.sorted.bam" 
	threads: maxthreads
	output:
		temp("sv/{sample}.sniffles.unsorted.vcf")
	log:
		"sv/{sample}.sniffles.log"
	shell: 
		"sniffles -m {input} -v {output} -t {threads} "
		"--tmp_file sv/{wildcards.sample}.sniffles.tmp " 
		"-s 2  &> {log}"

rule sort_vcf:
	input:
		"sv/{sample}.unsorted.vcf"
	output:
		"sv/{sample}.sorted.vcf"
	shell:
		"bcftools sort -o {output} -T sv {input}"

rule merge_sv_before_forcecall:
	input:
		expand("sv/{sample}.minimap2.sniffles.sorted.vcf",
			sample=nanopore_tb['sample']) # change "minimap2"/"ngmlr" to change aligner
	output:
		raw=temp("sv/raw_calls.txt"),
		vcf="sv/SURVIVOR_merged_1kbpdist_typesave.vcf"
	shell:
		"echo {input} | tr \" \" \"\\n\" > {output.raw} && "
		"SURVIVOR merge {output.raw} 1000 1 1 -1 -1 -1 {output.vcf}"

# https://github.com/fritzsedlazeck/Sniffles/wiki/SV-calling-for-a-population
rule sniffles_sv_forcecall:
	input:
		bam="bam/{sample}.sorted.bam",
		candidates="sv/SURVIVOR_merged_1kbpdist_typesave.vcf"
	threads: maxthreads
	output:
		temp("sv/{sample}.SURVIVOR.unsorted.vcf")
	log:
		"sv/{sample}.SURVIVORsniffles.log"
	shell:
		"sniffles -m {input.bam} -v {output} -t {threads} "
		"--Ivcf {input.candidates} --tmp_file "
		"sv/{wildcards.sample}.SURVIVORsniffles.tmp " 
		"-s 2 -n -1 --genotype --cluster &> {log}"

rule merge_sv_final:
	input:
		expand("sv/{sample}.minimap2.SURVIVOR.sorted.vcf",
			sample=nanopore_tb['sample']) # change "minimap2"/"ngmlr" to change aligner
	output:
		raw="sv/SURVIVORsniffles_calls.txt",
		vcf="sv/merged_final_SURVIVOR_1kbpdist_typesave.vcf"
	shell:
		"echo {input} | tr \" \" \"\\n\" > {output.raw} && "
		"SURVIVOR merge {output.raw} 1000 -1 1 -1 -1 -1 {output.vcf}"
