#!/usr/bin/Rscript
library(tidyverse)
library(GenomicRanges)
library(parallel)
library(getopt)
dir=dirname(get_Rscript_filename())
if (is.na(NA)) dir = "."
source(file.path(dir,"../../../methylation/methylation_plot_utils.R"))
cores = ceiling(detectCores()/2)

pdpath = "/shared/Data/chosigma_samples.txt"
mfreqdir = "/dilithium/Data/transgene"
regpath = "/dilithium/Data/NGS/projects/choSigmaAgingStudy/insert.bed"
plotdir = "/dilithium/Data/plots"

reg = load_db(regpath)

pd = read_tsv(pdpath)
pd$day = factor(pd$day)
pd$filepath = file.path(mfreqdir,paste0(pd$samples,".transgene.flank.meth.freq.txt"))

cnames = c("chrom","start","strand","meth","unmeth","dinuc","trinuc")
mfreq.list = mclapply(seq_along(pd$filepath),mc.cores=cores,function(i){
    read_tsv(pd$filepath[i],cnames)%>%
        mutate(sample=pd$samples[i])
})

mfreq = do.call(rbind,mfreq.list)
mfreq = bind_cols(mfreq,pd[match(mfreq$sample,pd$samples),-1]) %>%
    mutate(lab = paste(cell,glut)) %>%
    mutate(cov=meth+unmeth,freq=meth/cov) %>%
    filter(cov>1)

plotpath = file.path(plotdir,"methylation_transgene_flank_all.pdf")
pdf(plotpath,width=6,height=4,useDingbats=F)

g = ggplot(mfreq,aes(x=start,y=freq,color=lab,group=lab))+
    facet_wrap(.~day)+
    geom_point(size=0.5,alpha=0.5)+geom_smooth(size=0.5,se=F,span=0.7)+
    lims(x=c(15304000,15306000),y=c(0,1))+
    geom_vline(xintercept=15305200,linetype=2)+
    theme_bw()+
    labs(x="Coordinate along scaffold_16",y="Average methylation")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"),
          legend.position="bottom")
print(g)

dev.off()
