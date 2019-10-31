# AMBIC Epigenomics

This is the git repository for AMBIC Epigenomics project with Betenbaugh and Timp lab at JHU.  
Use this repository for processing the raw data prior to analysis.  
We are working on uploading scripts that are being used for our analysis of the data.

## Getting Started

### Prerequisites

The list of packages necessary for nanopore analysis :
+ [minimap2](https://github.com/lh3/minimap2)
+ [samtools](http://www.htslib.org/download/)
+ [htslib](http://www.htslib.org/)
+ [nanopolish](https://github.com/jts/nanopolish)
+ [bedtools](https://bedtools.readthedocs.io/en/latest/)
+ [wigToBigWig (from UCSC toolkit)](https://genome.ucsc.edu/util.html)
+ [sniffles](https://github.com/fritzsedlazeck/Sniffles)
+ [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
+ [bcftools](http://www.htslib.org/download/)
+ [pysam](https://pysam.readthedocs.io/en/latest/api.html)  

R packages : 
+ [bsseq](http://bioconductor.org/packages/release/bioc/html/bsseq.html)
+ [tidyverse](https://www.tidyverse.org/)

The list of packages necessary for ATAC-seq analysis :
+ [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
+ [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [samtools](http://www.htslib.org/doc/samtools.html)
+ [htslib](http://www.htslib.org/)
+ [picard](https://broadinstitute.github.io/picard/)
+ [bedtools](https://bedtools.readthedocs.io/en/latest/)
+ [macs2](https://github.com/taoliu/MACS)
+ [subread](http://bioinf.wehi.edu.au/subread-package/)

### Installing

All of these tools are available via conda (environment.yml) or via the website links.  
To install via conda, 
```
conda env update -n base --file ./environment.yml # to install to the base environment
conda env create -f ./environment.yml # to create a new environment
```

## Snakemake pipeline

We are using [snakemake](https://snakemake.readthedocs.io/en/stable/) to generate reproducible and easy-to-follow pipelines.  
You can either install snakemake and follow the steps to run the pipeline or directly look at the snakemake scripts to determine what commands to perform.  
Note : The snakemake workflow may not work completely as we have tested the workflow in limited settings. You may have to make edits to the code to make it work.

### Preparing snakemake

To use snakemake, first modify the **snakemake_config.yml** file as appropriate for your environment.

Secondly, modify the **nanopore_sample_info.csv** and **atacseq_sample_info.csv** as necessary.  

For nanopore sequencing data, we will start from basecalled data;  
set up the data as follows :

* fastq name : [sample].fastq.gz 
* summary name (optional) : [sample].summary.txt
* move fastq and summary files to [workdir]/reads
* fast5 in [workdir]/reads/[sample]/

For ATAC-seq data :

* fastq name : [sample].fastq.gz
* move fastq files into [workdir]/data/atacseq/fastq/

### Running snakemake

To run snakemake, use the snakemake command in the root directory of this repo.
```
snakemake --cores 16
```

For running specific parts of the pipeline, any of the following (currently broken): 
```
snakemake --cores 16 parse_nanopore
snakemake --cores 16 nanopore_methylation
snakemake --cores 16 nanopore_sv
snakemake --cores 16 parse_atacseq
```

For specific outputs,
```
snakemake --cores 16 path/to/output
```
e.g. :
```
snakemake --cores 16 methylation/[sample].cpg.meth.tsv.gz
```
will perform methylation calling on [sample] and any steps before that automatically, using 16 cores.