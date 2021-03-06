library(tidyverse)
library(reshape2)
plotdir="/home/isac/Dropbox/Data/ambic/atac/plots"

## insertions
plasmid="4_0cdhfr_vrc01wtg1m3_dgv"
idir="/dilithium/Data/Nanopore/Analysis/171025_cho/insertion/bed"
isuffix=".insertion.bed"
ipaths=system(paste0("ls ",idir,"/*",isuffix),intern=T)
##reps and labs will be used further down
labs=sapply(strsplit(sapply(strsplit(ipaths,"/"),"[[",9),"[.]"),"[[",1)
reps=sapply(strsplit(labs,"choNIH"),"[[",2)
## this pd should be used
pd=tibble(lab=reps,pheno=rep(c("host","IgG"),each=3),rep=rep(1:3,2),nanopore=labs)

cnames=c("chr","start","end","readname","score","strand")
inlist=lapply(ipaths,read_tsv,col_names=cnames)

bin=1000
win=1000
insertions=do.call(rbind,inlist)%>%filter(chr!=plasmid)
ins.bins=transmute(insertions,
                   chr=chr,
                   start=round(start/bin)*bin-win,
                   end=round(end/bin)*bin+win) %>%
    distinct(chr,start,end,.keep_all=T) %>%
    mutate_if(is.numeric,as.integer)
write_tsv(x=ins.bins,path=file.path(idir,"insertions.bed"))

## atac-seq data
covpath="/dilithium/Data/NGS/Aligned/171025_choatac/insertion/insertion.cov.txt"
ins.cov.all=read_tsv(covpath)
cov.tot=c(27.6,29.9,28.5,31.3,34.3,30.1)
normfactor=min(cov.tot)/cov.tot
ins.cov.all[,-c(1:3)]=round(ins.cov.all[,-c(1:3)]*normfactor)
ins.cov.all$id=seq(length(ins.cov.all$chr))
pd$atacseq=names(ins.cov.all)[-c(1:3)]

## data plotting
require(reshape2)
cov.melt=melt(ins.cov.all[,-c(1:3)],id.vars="id")
cov.melt$samp=rep(pd$pheno,each=dim(ins.cov.all)[1])
g=ggplot(data=cov.melt,mapping=aes(x=variable,y=value,color=samp,fill=samp,group=variable))+
    theme_bw()
g.viol=g+geom_violin(alpha=0.1,draw_quantiles=c(0.25,0.5,0.75))
g.den=g+geom_density(mapping=aes(x=value,y=..density..,
                                    color=samp,fill=samp,group=variable),
                        alpha=0.1)
g.box=g+geom_boxplot(alpha=0.1,outlier.shape=NA)
require(ggridges)
g.joy=g+geom_density_ridges(mapping=aes(x=value,y=variable),alpha=0.2,scale=6)
pdf(file.path(plotdir,"cov_density.pdf"),width=8,height=6)
print(g.viol)
print(g.den)
print(g.box)
print(g.joy)
dev.off()
##distance matrix
cov.vals=ins.cov.all[,-c(1:3)]
distance=dist(t(cov.vals))
clusters=hclust(distance)
cov.cor=cor(ins.cov.all[,-c(1:3)],method="pearson")

##function for only lookin at triangles
Lowertri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
Uppertri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

dist.mat = as.matrix(distance)
idx=order(rownames(dist.mat))
idx = match(labels(as.dendrogram(clusters)),rownames(dist.mat))
heat.mat = dist.mat[idx,idx]
rownames(heat.mat)=colnames(heat.mat)=seq(length(heat.mat[,1]))
heatmat.up=Uppertri(heat.mat)
dist.melt = melt(heatmat.up,na.rm=TRUE)

rownames(cov.cor)=colnames(cov.cor)=seq(length(heat.mat[,1]))
cor.low=round(Lowertri(cov.cor),3)
cor.melt=melt(cor.low,na.rm=TRUE)


q = ggplot(dist.melt,aes(x=Var1,y=Var2))+geom_tile(aes(fill=value))+
    geom_text(data=cor.melt,aes(x=Var1,Var2,label=value,size=10),color="black")+
        scale_fill_gradient(low="black",high="white",na.value="transparent")+
            theme_bw()+coord_fixed()+
                theme(legend.position="none",
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.ticks=element_blank(),
                      axis.title.x=element_blank(),
                      axis.text=element_blank(),
                      axis.title.y=element_blank())

pdf(file.path(plotdir,"atac_dendrogram.pdf"))
plot(clusters)
print(q)
dev.off()

##scatterplots
cov.spread=as.tibble(t(spread(cov.melt,key=id,value=value)[,-1]))
names(cov.spread)=pd$lab
pdf(file.path(plotdir,"covCorrelationPlot.pdf"),useDingbats=F,width=8,height=6)
require(GGally)
print(ggpairs(data=cov.spread,ggplot2::aes(alpha=0.1),
              lower=list(continuous = wrap("points",alpha=0.2,size=0.5)))+
      ggplot2::theme_bw())
dev.off()

