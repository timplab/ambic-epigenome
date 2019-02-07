#!/usr/bin/Rscript

library(tidyverse)
library(reshape2)


## data path
fin="/atium/Data/NGS/Aligned/170718_ambicatac/manual/insertion_coverage.txt"
wdir=dirname(fin)
fname=basename(fin)
base=strsplit(fname,"[.]")[[1]][1]

## read in data and get length
cnames=c("seqname","start","end","name","cov1","cov2","cov3")
cov=read_tsv(fin,col_names=cnames) %>%
    mutate(cov=cov1+cov2+cov3,
           relcov=cov/cov[5])
cov=cov[3:5,]
gcov=ggplot(data=cov,mapping=aes(x=name,y=relcov,fill=name))+theme_bw()+
    geom_col()+#    scale_y_log10()+
    guides(fill=FALSE)+
    labs(title="ATAC-seq coverage",
         x="Region")

outname=file.path(wdir,paste0(base,"_insertionCov.pdf"))
pdf(outname,width=3,height=2)
print(gcov)
dev.off()
