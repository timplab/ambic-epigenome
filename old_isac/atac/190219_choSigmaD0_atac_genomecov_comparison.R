#!/usr/env Rscript
library(tidyverse)
dir="/kyber/Data/NGS/projects/ambic/agingstudy/atacseq/190211_choSigmaD0_atac/data/atacseq/genomecov"
outdir="~/Dropbox/tmp"
fp = system(paste("find",dir,"-name \"*genomecov.bed\"",sep=" "),intern=T)
fn = sapply(strsplit(fp,"/"),"[[",13)
base = sapply(strsplit(fn,"[.]"),"[[",1)
pd = tibble(base=base,fp=fp)

controlpath = "/dilithium/Data/NGS/projects/gm12878/bsseq/bed/bsseq.genomecov.bed"

dat.list = lapply(seq_along(pd$base),function(i){
    read_tsv(pd$fp[i],col_names=F) %>%
        mutate(sample=pd$base[i]) %>%
        filter(X1=="genome") %>%
        rename(cov = X2,frac = X5) %>%
        select(cov,frac,sample)
})

control = read_tsv(controlpath,col_names=F) %>%
    rename(cov=X2,frac=X5) %>%
    filter(X1=="genome") %>%
    select(cov,frac)
controlcov = sum(control$cov * control$frac)

control.list = lapply(dat.list,function(x){
    totcov = sum(x$cov * x$frac)
    bin = controlcov/totcov
    control %>%
        mutate(cov = round(cov/bin)) %>%
        group_by(cov)%>%summarize(frac = sum(frac)) %>%
        mutate(sample="control")
})

outpath=file.path(outdir,"genomecov.pdf")
pdf(outpath,useDingbats=F)
#g = ggplot(control,aes(x=cov,y=frac))+
#    geom_line()+
#    scale_y_log10()
#print(g)
for (i in seq_along(dat.list)){
    dat.plt = dat.list[[i]]
##    control.sub = control.list[[i]]
##    dat.plt = do.call(rbind,list(dat.sub,control.sub))
    g = ggplot(dat.plt,aes(x=cov,y=frac,group=sample,color=sample))+
        geom_line()+
        theme_bw()
    g.log = g + scale_y_log10()
    g.high = g + lims(y=c(0,1e-6))
    print(g.log)
    print(g.high)
}
dev.off()
