# AMBIC Epigenomics

This is the git repository for AMBIC Epigenomics project with Betenbaugh and Timp lab at JHU.

This repo is still being put together, so it will only contain parts of the analysis.

At the current stage, this should provide enough information for nanopore methylation data parsing

## Getting Started

We are implementing the pipeline using the snakemake package.

### Prerequisites

The list of packages necessary for nanopore analysis :
```
snakemake
ngmlr
samtools
htslib
nanopolish
bedtools
wigToBigWig (from UCSC toolkit)
```

The list of packages necessary for ATAC-seq analysis :
```
trim_galore
bowtie2
samtools
htslib
picard
bedtools
macs2
subread
```

### Installing

## Analysis pipeline

### Preparing snakemake

To use snakemake, first modify the snakemake_config.yml file as appropriate for your environment

Secondly, modify the nanopore_sample_info.csv as necessary

Then move data as follows (relative to "workdir") :

* fastq and summary files in data/nanopore/reads
* fast5 tarball in data/nanopore/tar or untarred into data/nanopore/reads/[sample_name]

### Running snakemake

To run snakemake, use the snakemake command in the root directory of the repo.
```
snakemake --cores 16
```

For running specific parts of the pipeline,
```
snakemake --cores 16 parse_nanopore
snakemake --cores 16 parse_atacseq
```
