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
mfreqdir = "/dilithium/Data/mfreq_pooled"
regpath = "/dilithium/Data/NGS/projects/choSigmaAgingStudy/insert.bed"
plotdir = "/dilithium/Data/plots"

reg = load_db(regpath)

pd = read_tsv(pdpath)
pd$day = factor(pd$day)
pd$filepath = file.path(mfreqdir,paste0(pd$samples,".chok1gshd_sigmaIgG.methfreq.txt.gz"))

mfreq.list = mclapply(seq_along(pd$filepath),mc.cores=cores,function(i){
    tabix_mfreq(pd$filepath[i],reg,cov=1,trinuc_exclude=NULL)%>%
        mutate(sample=pd$samples[i])
})

mfreq = do.call(rbind,mfreq.list)
mfreq = bind_cols(mfreq,pd[match(mfreq$sample,pd$samples),-1]) %>%
    mutate(lab = paste(cell,glut))
mfreq = bind_rows(mfreq,
                  x = tibble(chrom = c("ambic_sigma_IgG_HC","ambic_sigma_IgG_LC"),
                             start = 0,strand = "*",meth=0,
                             unmeth=10,dinuc="CG",trinuc="CCG",cov=10,freq=0,
                             end=0,sample="choSigmaHostD0",cell="Host",day=factor(0),
                             glut="Glut",filepath="place",lab="Host Glut")) # faux data for coloring

plotpath = file.path(plotdir,"methylation_transgene_body_all.pdf")
pdf(plotpath,width=6,height=4,useDingbats=F)
for (chrom in unique(mfreq$chrom)){
    mfreq.sub = mfreq[which(mfreq$chrom==chrom),]
    g = ggplot(mfreq.sub,aes(x=start,y=freq,color=lab,group=lab))+
        facet_wrap(.~day)+
        geom_point(size=0.5,alpha=0.5)+geom_smooth(size=0.5,se=F,span=0.7)+
        theme_bw()+
        labs(x=paste0("Coordinate along ",chrom),y="Average methylation")+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(color="black"),
              legend.position="bottom")
    print(g)
}
dev.off()
