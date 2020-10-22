#!/usr/bin/Rscript

library(tidyverse)
library(GenomicRanges)

root="/dilithium/Data/Nanopore/Analysis/171025_cho"
wdir="/dilithium/Data/Nanopore/Analysis/171025_cho/meth"
atacdir="/dilithium/Data/NGS/Aligned/171025_choatac/insertion"
outdir="~/Dropbox/Data/ambic/insertion"

insertion.fh=file.path(root,"cho.insertion.bed")
host.atac.fh=file.path(atacdir,"host.insert.1kb.bed")
igg.atac.fh=file.path(atacdir,"IgG.insert.1kb.bed")
host.fh=file.path(wdir,"choNIHhost.insertion.methfreq.tsv")
igg.fh=file.path(wdir,"choNIHIgG.insertion.methfreq.tsv")


insertion=read_tsv(file=insertion.fh,col_names=c("chr","start","end")) %>%
    mutate(start=end-2000)
host=read_tsv(file=host.fh)
igg=read_tsv(file=igg.fh)
bedcov.cnames=c("chr","start","end","count","bpcover","bptotal","coverfreq")
host.atac=read_tsv(file=host.atac.fh,col_names=bedcov.cnames)
igg.atac=read_tsv(file=igg.atac.fh,col_names=bedcov.cnames)
ins.gr=GRanges(insertion)
host.gr=GRanges(host)
igg.gr=GRanges(igg)

host.sum=findOverlaps(ins.gr,host.gr) %>%
    as.data.frame() %>%
    mutate(meth=host.gr$methylated_frequency[subjectHits]) %>%
    group_by(queryHits) %>%
    summarize(mean=mean(meth),
              num=n())

igg.sum=findOverlaps(ins.gr,igg.gr) %>%
    as.data.frame() %>%
    mutate(meth=igg.gr$methylated_frequency[subjectHits]) %>%
    group_by(queryHits) %>%
    summarize(mean=mean(meth),
              num=n())

meth.sum=insertion %>%
    mutate(hostmeth=host.sum$mean,
           iggmeth=igg.sum$mean,
           hostatac=host.atac$count,
           iggatac=igg.atac$count)

write_tsv(x=meth.sum,path=file.path(outdir,"insertion.sum.tsv"))
