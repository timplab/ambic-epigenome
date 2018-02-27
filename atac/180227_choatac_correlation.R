#!/usr/bin/Rscript

library(tidyverse)
plotdir="~/Dropbox/Data/Atac-seq/180223_choOmniAtac/plot"
omnidir="/dilithium/Data/NGS/Aligned/180223_choOmniAtac/tss"
atacdir="/dilithium/Data/NGS/Aligned/171025_choatac/bed/tss"
omnipaths=system(intern=T,command=paste0('find ',omnidir," -name *txt"," -type f"))
atacpaths=system(intern=T,command=paste0('find ',atacdir," -name *txt"," -type f"))
omnilabs=sapply(strsplit(sapply(strsplit(sapply(strsplit(omnipaths,"/"),"[[",8),"_"),"[[",1),"[.]"),"[[",1)


cnames=c("distance","gene")
omni.list <- lapply(omnipaths,function(x)read_tsv(x,col_names=cnames))
atac.list <- lapply(atacpaths,function(x)read_tsv(x,col_names=cnames))
atac=atac.list[[4]]
b=10
summarizeCov <- function(x,b){
    x%>%mutate(bindist=round(distance/b)*b)%>%
        group_by(bindist,gene)%>%
        summarize(cov=n(),
                  covnorm=cov/dim(x)[1]*100)
}

omni.sum <- omni.list%>%
    map(function(x){summarizeCov(x,b)})

atac.sum <- summarizeCov(atac,b) 

cov.comp <- omni.sum %>%
    map(function(x){
        list(x,atac.sum)%>%
            reduce(inner_join,by=c("bindist","gene"))%>%
            group_by(bindist,gene)%>%
            select(starts_with("covnorm"))%>%
            filter(covnorm.x<1)
    })

comp.cor <- sapply(cov.comp,function(x){
    cor(x$covnorm.x,x$covnorm.y)})
comp.lm <- lapply(cov.comp,function(x){
    lm(covnorm.x~covnorm.y,x)})
glist=lapply(cov.comp,function(x){
    ggplot(x,aes(x=covnorm.x,y=covnorm.y))+
        geom_point(size=0.5)+theme_bw()+
        labs(x="Omni-ATAC",y="ATAC-seq",title="Coverage Scatterplot")
})

pdf(file.path(plotdir,"covCorrelation.pdf"),useDingbats=F)
for (i in seq_along(glist)){
    print(glist[[i]]+ggtitle(omnilabs[i]))
}
dev.off()
