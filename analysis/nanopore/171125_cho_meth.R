library(tidyverse)
library(reshape2)
fdir="/dilithium/Data/Nanopore/Analysis/171025_cho/meth"
suffix=".methfreq.tsv"
paths=system(paste0("ls ",fdir,"/*",suffix),intern=T)
plotdir="/home/isac/Dropbox/Data/ambic/nanopore/plots"

labs=sapply(strsplit(sapply(strsplit(paths,"/"),"[[",8),"[.]"),"[[",1)
meths=lapply(paths,read_tsv)
names(meths)=labs

makeBin=function(x,bin=10000){
    mutate(x,
             start=round(start/bin)*bin,
             end=round(end/bin)*bin) %>%
        group_by(chromosome,start)%>%
        summarize(calls=sum(called_sites),
                  methylated=sum(called_sites_methylated))%>%
        mutate(freq=methylated/calls)
}
addLabCols=function(l){
    lapply(seq_along(l),function(i){
        l[[i]]%>%mutate(lab=names(l)[i])})}

#x = makeBin(meths[[1]])
b=100
bins=lapply(meths,function(x) makeBin(x,b))
bins=addLabCols(bins)

##find overlaps
binovl=Reduce(function(df1,df2)inner_join(df1,df2,by=c("chromosome","start")),bins)
cnames=colnames(binovl)
indfreq=grep("freq",cnames)

freqmat=binovl[,indfreq]
colnames(freqmat)=labs

##let's make violin or density plot of methylation:
len=function(x){dim(x)[1]}
meth.melt=melt(freqmat,id.vars=NULL,variable.name="id",value.name="Methylation")
samp=rep(c("host","IgG"),each=3)
meth.melt$sample=rep(samp,each=len(freqmat))
g=ggplot(data=meth.melt,
         mapping=aes(x=id,y=Methylation,color=sample,fill=sample,group=id))+
    theme_bw()
g.viol=g+geom_violin(alpha=0.1,draw_quantiles=c(0.25,0.5,0.75))
g.den=g+geom_density(mapping=aes(x=Methylation,y=..density..,
                                 color=sample,fill=sample,group=id),
                     alpha=0.1)
g.box=g+geom_boxplot(alpha=0.1,outlier.shape=NA)
pdf(file.path(plotdir,"meth_density.pdf"),width=8,height=6)
print(g.viol)
print(g.den)
print(g.box)
dev.off()

##distance matrix?
distance = dist(t(freqmat))
clusters = hclust(distance)
meth.cor=cor(freqmat,method="pearson")

##function for only lookin at triangles
# Get lower triangle of the correlation matrix
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

##pearson correlation
rownames(meth.cor)=colnames(meth.cor)=seq(length(heat.mat[,1]))
methcor.low=round(Lowertri(meth.cor),2)
cor.melt=melt(methcor.low,na.rm=TRUE)


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

pdf(file.path(plotdir,"nanopore_dendrogram.pdf"))
plot(clusters)
print(q)
dev.off()

