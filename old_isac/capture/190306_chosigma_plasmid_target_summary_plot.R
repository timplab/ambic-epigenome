#!/usr/bin/env Rscript
library(tidyverse)
library(wesanderson)

dir = "/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target"
outdir = "~/Dropbox/Data/ambic/capture/plots"
fapaths = system(paste("find",dir,"-name \"*summary.txt\""),intern=T)
samples = sapply(strsplit(basename(fapaths),"[.]"),"[[",1)

pd = tibble(sample=samples,fp=fapaths)

cnames = c("qname","qstart","qend","rname","rstart","rend","strand","qlen")
dat.list = lapply(seq_along(pd$sample),function(i){
    read_tsv(pd$fp[i],col_names=cnames)%>%
        mutate(sample=pd$sample[i])
})

# get the contig with most alignments
long.list = lapply(dat.list,function(x){
    contigs = names(sort(table(x$qname),T))[1:2]
    x[x$qname %in% contigs,]
})


# plot
dat = do.call(rbind,long.list) %>%
    mutate(label = paste0(rname,":",rstart,"-",rend))
yfac = factor(dat$qname)
levels(yfac) = c(1,3)
ymax = as.numeric(as.character(yfac))
cols = factor(dat$rname)

dat = dat %>%
    mutate(ymax=ymax,
           ymin=ymax-1,
           col=cols)

dat.contigs = dat %>%
    select(qname,qlen,sample,ymin,ymax) %>%
    unique()
    
g = ggplot(dat,aes(ymin=ymin,ymax=ymax,
                   xmin=qstart,xmax=qend,fill=rname))+
    geom_rect(inherit.aes=F,
              data=dat.contigs,
              mapping=aes(xmin=0,xmax=qlen,ymin=ymin,ymax=ymax))+
    geom_label(inherit.aes=F,nudge_y = 0.2,
              data=dat.contigs,hjust = 0,
              mapping=aes(x=0,y=ymax,label=sample))+
    geom_rect()+
    scale_fill_manual(values = wes_palette("Rushmore1"))+
    theme_bw() + lims(y=c(0,4))+
    theme(legend.position="bottom")

plotpath = file.path(outdir,"contig_summary.pdf")
pdf(plotpath,useDingbats=F,width=8,heigh=4)
print(g)
dev.off()
