#!/usr/bin/snakemake --snakefile
import pandas as pd
from snakemake.utils import validate
nanopore_tb = pd.read_csv(config['codedir']+"/nanopore_sample_info.csv")

##################################################
#
# alignment
#
##################################################

rule ngmlr_align:
	input:
		"{dir}/reads/{sample}.fastq.gz"
	params:
		config['reference']
	threads: maxthreads
	output:
		"{dir}/bam/{sample}.ngmlr.sorted.bam"
	log:
		"{dir}/bam/{sample}.ngmlr.align.log"
	shell:
		"ngmlr -x ont --bam-fix "
		"-t {threads} -r {params} -q {input} 2> {log} | "
		"samtools view -q 20 -b - | "
		"samtools sort -T {wildcards.dir}/bam/{wildcards.sample}.sorting "
		"-o {output} && "
		"samtools index {output}"
	
rule minimap2_makeindex:
	input:
		config['reference']
	output:
		os.path.splitext(config['reference'])[0] + ".mapont.mmi"
	shell:
		"minimap2 -x map-ont -d {output} {input}"

rule minimap2_align:
	input:
		os.path.splitext(config['reference'])[0] + ".mapont.mmi",
		"{dir}/reads/{sample}.fastq.gz"
	threads: maxthreads
	output:
		"{dir}/bam/{sample}.minimap2.sorted.bam"
	log:
		"{dir}/bam/{sample}.minimap2.align.log"
	shell:
		"minimap2 --MD -L -t {threads} -a {input} 2> {log} | "
		"samtools view -q 20 -b - | "
		"samtools sort -T {wildcards.dir}/bam/{wildcards.sample}.sorting "
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
		fq="{dir}/reads/{sample}.fastq.gz",
		bam="{dir}/bam/{sample}.ngmlr.sorted.bam", # change "minimap2"/"ngmlr" to change aligner
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

##################################################
#
# SV detection 
#
##################################################

rule sniffles_sv:
	input:
		"{dir}/bam/{sample}.sorted.bam" 
	threads: maxthreads
	output:
		temp("{dir}/sv/{sample}.sniffles.unsorted.vcf")
	log:
		"{dir}/sv/{sample}.sniffles.log"
	shell: 
		"sniffles -m {input} -v {output} -t {threads} "
		"--tmp_file {wildcards.dir}/sv/{wildcards.sample}.sniffles.tmp " 
		"-s 5 -n -1 --genotype --cluster &> {log}"

rule sort_vcf:
	input:
		"{dir}/sv/{sample}.unsorted.vcf"
	output:
		"{dir}/sv/{sample}.sorted.vcf"
	shell:
		"awk '/#/ {{print;next}}{{ print | \"sort -k1,1 -k2,2n \" }}' "
		"{input} > {output}"

rule merge_sv_before_forcecall:
	input:
		expand("data/nanopore/sv/{sample}.minimap2.sniffles.sorted.vcf",
			sample=nanopore_tb['sample']) # change "minimap2"/"ngmlr" to change aligner
	output:
		raw=temp("data/nanopore/sv/raw_calls.txt"),
		vcf=temp("data/nanopore/sv/SURVIVOR_merged_1kbpdist_typesave.vcf")
	shell:
		"echo {input} | tr \" \" \"\\n\" > {output.raw} && "
		"SURVIVOR merge {output.raw} 1000 1 1 -1 -1 -1 {output.vcf}"

rule sniffles_sv_forcecall:
	input:
		bam="{dir}/bam/{sample}.sorted.bam",
		candidates="data/nanopore/sv/SURVIVOR_merged_1kbpdist_typesave.vcf"
	threads: maxthreads
	output:
		temp("{dir}/sv/{sample}.SURVIVOR.unsorted.vcf")
	log:
		"{dir}/sv/{sample}.SURVIVORsniffles.log"
	shell:
		"sniffles -m {input.bam} -v {output} -t {threads} "
		"--Ivcf {input.candidates} --tmp_file "
		"{wildcards.dir}/sv/{wildcards.sample}.SURVIVORsniffles.tmp " 
		"-s 5 -n -1 --genotype --cluster &> {log}"

rule merge_sv_final:
	input:
		expand("data/nanopore/sv/{sample}.minimap2.SURVIVOR.sorted.vcf",
			sample=nanopore_tb['sample']) # change "minimap2"/"ngmlr" to change aligner
	output:
		raw="data/nanopore/sv/SURVIVORsniffles_calls.txt",
		vcf="data/nanopore/sv/merged_final_SURVIVOR_1kbpdist_typesave.vcf"
	shell:
		"echo {input} | tr \" \" \"\\n\" > {output.raw} && "
		"SURVIVOR merge {output.raw} 1000 -1 1 -1 -1 -1 {output.vcf}"
