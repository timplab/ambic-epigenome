#!/usr/bin/snakemake --snakefile
# referenced : https://github.com/ENCODE-DCC/atac-seq-pipeline
import pandas as pd
from snakemake.utils import validate 
atacseq_tb = pd.read_csv(config['codedir']+"/atacseq_sample_info.csv")

rule trim_fastq: 
	input: 
		"{dir}/fastq/{sample}.fastq.gz" 
	output: 
		"{dir}/fastq_trimmed/{sample}_trimmed.fq.gz"
	log: 
		"{dir}/log/trim_fastq/{sample}.trim.log" 
	shell: 
		"trim_galore {input} "
		"-o {wildcards.dir}/fastq_trimmed &> {log}"

rule index_bowtie2:
  input:
    config['reference']
  params:
    config['reference'].split(".fa")[0]
  threads:
    maxthreads
  output:
    config['reference'].split(".fa")[0] + ".1.bt2"
  log:
    config['reference'].split(".fa")[0] + ".btidx.log"
  shell:
    "bowtie2-build {input} {params} --threads {threads} &> {log}"

rule atacseq_align: 
	input: 
		fq="{dir}/fastq_trimmed/{sample}_trimmed.fq.gz",
		btidx=config['reference'].split(".fa")[0] + ".1.bt2"
	params: 
		config['reference'].split(".fa")[0] 
	threads: 
		maxthreads 
	output: 
		"{dir}/bam/{sample}.bowtiealign.bam"
	log: 
		"{dir}/log/align/{sample}.align.log" 
	shell: 
		"bowtie2 -k 4 -p {threads} -t --local " 
		"-x {params} -U {input.fq} 2> {log} | " 
		"samtools view -Sb - > {output}"

rule filter_multiple_alignments: 
	input: 
		"{dir}/bam/{sample}.bowtiealign.bam"
	params:
		config['codedir']
	output:
		"{dir}/bam/{sample}.multimapfiltered.bam"
	shell:
		"samtools sort -n -T "
		"{wildcards.dir}/bam/{wildcards.sample}.qsorting {input} | "
		"samtools view -h - | "
		"python {params}/scripts/assign_multimappers.py -k 4 | "
		"samtools view -F 1804 -Su - | "
		"samtools sort -T "
		"{wildcards.dir}/bam/{wildcards.sample}.mpsorting "
		"-o {output}"

rule mark_duplicates:
	input:
		"{dir}/bam/{sample}.multimapfiltered.bam"
	output:
		temp("{dir}/bam/{sample}.dupmark.bam")
	threads:
		maxthreads
	log:
		"{dir}/log/removedup/{sample}.dupremoval.log"
	shell:
		"picard -Xmx8G -Xms256M -XX:ParallelGCThreads={threads} "
		"-Djava.io.tmpdir={wildcards.dir}/tmp "
		"MarkDuplicates INPUT={input} OUTPUT={output} "
		"METRICS_FILE={wildcards.dir}/bam/{wildcards.sample}.dupmetrics.txt "
		"VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true "
		"REMOVE_DUPLICATES=false &> {log}"

rule remove_duplicates:
	input:
		"{dir}/bam/{sample}.dupmark.bam"
	output:
		"{dir}/bam/{sample}.nodup.sorted.bam"
	shell:
		"samtools view -F 1804 -b {input} > {output} && "
		"samtools index {output}"

rule atacseq_bam_to_bed:
	input:
		"{dir}/bam/{sample}.nodup.sorted.bam"
	params:
		"{dir}/tmp"
	output:
		"{dir}/bed/{sample}.bed.gz"
	shell:
		"[ -e {params} ]||mkdir {params} && "
		"bedtools bamtobed -i {input} | "
		"awk 'BEGIN{{OFS=\"\t\"}}{{$4=\"N\";$5=\"1000\";print $0 }}' | "
		"sort -k1,1 -k2,2n -T {params} | "
		"bgzip > {output}"

rule atacseq_merge_bed:
	input:
		expand("{dir}/bed/{sample}.bed.gz",
			dir="{dir}",sample=atacseq_tb['sample'])
	params:
		"{dir}/tmp"
	output:
		"{dir}/bed/allsamples.bed.gz"
	shell:
		"[ -e {params} ]||mkdir {params} && "
		"gunzip -c {input} | "
		"sort -k1,1 -k2,2n -T {params} | "
		"bgzip > {output}"

rule find_peaks_from_pool:
	input:
		"{dir}/bed/allsamples.bed.gz"
	params:
		config['codedir']
	output:
		"{dir}/peaks/allsamples.peaks.saf"
	log:
		"{dir}/log/allsamples.peaks.log"
	shell:
		"macs2 callpeak -t {input} -f BED "
		"-n {wildcards.dir}/peaks/allsamples "
		"-g mm -p 0.01 --shift -75 --extsize 150 "
		"--nomodel -B --SPMR --keep-dup all --call-summits &> {log} && "
		"{params}/scripts/atacseq_bed_to_saf.sh "
		"{wildcards.dir}/peaks/allsamples_peaks.narrowPeak > {output}"

rule get_feature_counts:
	input:
		reads=expand("{dir}/bam/{sample}.nodup.sorted.bam",
			dir="{dir}",sample=atacseq_tb['sample']),
		peaks="{dir}/peaks/allsamples.peaks.saf"
	threads:
		maxthreads
	output:
		"{dir}/peaks/individual_peaks_counts.txt"
	log:
		"{dir}/log/individual.featurecounts.log"
	shell:
		"featureCounts --verbose -a {input.peaks} "
		"-o {output} -F SAF -T {threads} {input.reads} &> {log}"
