library(tidyverse)
library(reshape2)
plotdir="/home/isac/Dropbox/Data/ambic/integrative-analysis/plots"

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
pd$atacseq=names(ins.cov.all)[-c(1:3)]

## methylation data
mdir="/dilithium/Data/Nanopore/Analysis/171025_cho/meth"
methsuffix=".methfreq.tsv"
methpaths=sapply(labs,function(x) paste0(mdir,"/",x,methsuffix))
meths=lapply(methpaths,read_tsv)
names(meths)=labs

## bin methylation data
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

##make meth matrix
binovl=Reduce(function(df1,df2)inner_join(df1,df2,by=c("chromosome","start")),bins)
cnames=colnames(binovl)
indfreq=grep("freq",cnames)

freqmat=binovl[,c(1,2,indfreq)]
colnames(freqmat)=c("chr","pos",labs)

##find overlaps
require(GenomicRanges)
ins.gr=GRanges(ins.bins)
freq.gr=GRanges(seqnames=freqmat$chr,ranges=IRanges(start=freqmat$pos,width=1))
ins.ovl=findOverlaps(ins.gr,freq.gr)
ins.meth.all=as.tibble(t(sapply(seq_along(ins.bins$chr),function(i){
    ind=which(queryHits(ins.ovl)==i)
    meth=freqmat[subjectHits(ins.ovl)[ind],-c(1,2)]
    colMeans(meth)
})))

##get rid of NaNs
dropidx=which(is.na(ins.meth.all[,1]))
ins.meth=ins.meth.all[-dropidx,]
ins.cov=ins.cov.all[-dropidx,-c(1:3)]
ins.pos=ins.bins[-dropidx,]

##tidy the data for plotting
ins.cov$type="atac"
ins.meth$type="meth"
names(ins.cov)=names(ins.meth)=c(pd$lab,"type")
ins.cov$id=ins.meth$id=seq_along(ins.cov$type)
dat.df=rbind(ins.cov,ins.meth)
require(reshape2)
dat.melt=melt(dat.df,id.vars=c("type","id"))
dat.spread=spread(dat.melt,type,value)
dat.spread$samp=rep(pd$pheno)
dat.spread$rep=rep(pd$rep)
dat.sum=group_by(dat.spread,id,samp)%>%
    summarize(atac=sum(atac),
              meth=mean(meth)) %>%
    group_by(id)%>%summarize(atac=diff(atac)/sum(atac)*2,meth=diff(meth))
#dat.sum=group_by(dat.spread,id,rep) %>%
#    summarize(atac=diff(atac)/sum(atac)*2,meth=diff(meth))

g=ggplot(data=dat.spread,mapping=aes(y=meth,x=atac,color=samp))+theme_bw()
g.point=g + geom_point()
gdiff=ggplot(data=dat.sum,mapping=aes(y=meth,x=atac))+theme_bw()+
    geom_point()
g.atac=ggplot(data=dat.sum,mapping=aes(x=atac,y=..density..))+
    geom_density()+theme_bw()
g.meth=ggplot(data=dat.sum,mapping=aes(x=meth,y=..density..))+
    geom_density()+theme_bw()
require(ggridges)
g.atac.joy=g+geom_density_ridges(mapping=aes(x=atac,y=variable,fill=samp),alpha=0.2,scale=6)
g.meth.joy=g+geom_density_ridges(mapping=aes(x=meth,y=variable,fill=samp),alpha=0.2,scale=3)

pdf(file.path(plotdir,"insertion_difference_density.pdf"),width=8,height=6)
print(g.atac)
print(g.meth)
print(g.atac.joy)
print(g.meth.joy)
dev.off()

pdf(file.path(plotdir,"scatterplot.pdf"),width=8,height=6,useDingbats=F)
print(g.point)
print(gdiff)
dev.off()
