#!/usr/bin/env Rscript
library(tidyverse)
root = "/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target"
outdir = "~/Dropbox/Data/ambic/capture/plots"
fps = system(paste("find ",root,"-name \"*blast.txt\""),intern=T)
samples = sapply(strsplit(basename(fps),"_blast"),"[[",1)
pd = tibble(sample = samples, fp = fps)

blast_cnames = c("qname","rname","match","alen","mismatch","gap",
                 "qstart","qend","sstart","send","eval","bscore")
dat.list = lapply(seq_along(pd$fp),function(i){
    read_tsv(pd$fp[i],blast_cnames) %>%
        mutate(sample = pd$sample[i])
})

dat = do.call(rbind,dat.list)




dat.mostmatch = lapply(dat.list,function(x){
    cnt = as.tibble(table(x$qname[grepl("sigma",x$rname)]))
    read = cnt$Var1[which.max(cnt$n)]
    x[which(x$qname==read),]
})



