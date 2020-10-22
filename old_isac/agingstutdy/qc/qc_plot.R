#!/usr/bin/env Rscript
library(tidyverse)
root = "/kyber/Data/Nanopore/projects/ambic/sigma/qc"
qcfp = file.path(root,"190308_ambic_agingstudy_fastq_summary.csv")
outdir = "~/Dropbox/Data/ambic/aging_study/qc"

dat = read_csv(qcfp) %>%
    arrange(yieldGb) %>%
    filter(yieldGb>=10) %>%
    mutate(avglen=yieldGb*1e9/numreads)
