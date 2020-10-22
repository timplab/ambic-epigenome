#!/usr/bin/Rscript
library(tidyverse)

csvpath="/home/isac/Dropbox/Data/ambic/aging_study/Passage.csv"
outpath="/home/isac/Dropbox/Data/ambic/aging_study/plots/growth.pdf"
csv=read_csv(csvpath)

dat=csv%>%
    gather(Sample,PDL,-Passage)

g=ggplot(dat,aes(x=Passage,y=PDL,group=Sample,color=Sample))+
    geom_line()+geom_point(size=0.5)+
    theme_bw()+labs(y="Cumulative PDL",x="Passage Number")

pdf(outpath,width=6,height=3,useDingbats=FALSE)
print(g)
dev.off()

