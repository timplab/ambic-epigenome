#!/bin/bash

# using fastq generated from gilfunk
# He took out reads that have >= 400 bp alignment to the HC of IgG sigma transgene

# first get lengths

dir=/home/isac/Dropbox/Data/ambic/capture/analysis/
gunzip -c $dir/bothSIDES_363.fastq.gz | awk 'NR%4==0{ print length($1) }'  > $dir/bothSIDES_lengths.txt

# plot
knit_rmd.sh -d ~/Dropbox/Data/ambic/capture/analysis/ 191206_capture_plasmid_fastq_lengths.Rmd
