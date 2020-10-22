library(tidyverse)
beddir="/dilithium/Data/Nanopore/Analysis/171025_cho/insertion/bed"
bednames=system(paste0("ls ",beddir),T)
samps=sapply(strsplit(bednames,"[.]"),"[[",1)

pname="4_0cdhfr_vrc01wtg1m3_dgv"

cnames=c("chr","start","end","rname","score","strand")
beds=lapply(X=samps,FUN=function(x) {
    cbind(read_tsv(file=file.path(beddir,paste0(x,".insertion.bed")),col_names=cnames),x)})

beds.df=do.call(rbind.data.frame,beds)

beds.sum=beds.df %>%
    filter(chr!=pname) %>%
    mutate(len=end-start,bin=round(start/1000)) %>%
    filter(len>100) %>%
    group_by(chr,x,bin) %>%
    summarize(n()) %>%
    filter(grepl("IgG",x))

