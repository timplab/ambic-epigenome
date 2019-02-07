#!/usr/bin/Rscript
library(tidyverse)
library(parallel)
cores = ceiling(detectCores()/2)

root = "/dilithium/Data/NGS/projects/choSigmaAgingStudy"
pdpath = file.path(root,"atacseq_samples.txt")
covpath = file.path(root,"peaks/choSigmaAgingStudy_rep1_peaks_counts.txt")
plotdir = "/dilithium/Data/plots"

pd = read_tsv(pdpath)

cov = read_tsv(covpath,skip=1)
newcnames = sapply(strsplit(
        basename(names(cov)[7:length(names(cov))]),"[.]"),"[[",1)
names(cov)[7:length(names(cov))] = newcnames

# normalize
for (x in pd$samples) {
    cov.samp = cov[[x]]
    totcov = sum(cov.samp)
    cov.norm = cov.samp/totcov
    cov = mutate(cov,!!x := cov.norm)
}

# subset igg
cov.trans = cov[grep("ambic",cov$Chr),]
cov.gather = cov.trans %>%
    select(-Start,-End,-Strand,-Length,-Geneid) %>%
    gather("sample","coverage",-Chr)
cov.spread = cov.gather %>%
    spread(Chr,coverage)

cov.tb = bind_cols(pd[match(cov.spread$sample,pd$samples),],cov.spread)%>%
    mutate(totcov = ambic_sigma_IgG_LC+ambic_sigma_IgG_HC)

# plot
g = ggplot(cov.tb,aes(x=day,y=totcov,colour=paste(cell,glut),group=paste(cell,glut)))+
    geom_point()+geom_line(size=0.5)+theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color="black"))+
    labs(x="Day",y="Transgene coverage fraction")

pdf(file.path(plotdir,"181203_transgene_accessibility.pdf"),useDingbats=F,width=5,height=3)
print(g)
dev.off()
