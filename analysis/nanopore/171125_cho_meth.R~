library(tidyverse)
library(VennDiagram)

fdir="/dilithium/Data/Nanopore/Analysis/171025_cho/sniffles"
suffix=".sniffles.stringent.bedpe"
hostpath=paste0(fdir,"/choNIHhost",suffix)
iggpath=paste0(fdir,"/choNIHIgG",suffix)
plotdir="/home/isac/Dropbox/Data/ambic/nanopore/plots"

cnames=c("chrom1","start1","stop1","chrom2","start2","stop2","ID","score","strand1","strand2","type","number_of_split_reads","best_chr1","best_start","best_chr2","best_stop","len")
host=read_tsv(hostpath,col_names=cnames,skip=1)
igg=read_tsv(iggpath,col_names=cnames,skip=1)

makeBin=function(x,bin){
    y=mutate(x,
             start1=round(start1/bin)*bin,
             stop1=round(stop1/bin)*bin,
             start2=round(start2/bin)*bin,
             stop2=round(stop2/bin)*bin) %>%
        transmute(name=paste(chrom1,start1,stop1,chrom2,start2,stop2,strand1,strand2,type,sep="_"))%>%
        unlist(use.names=F)
    
}
                   
bin=1000
host.bin=makeBin(host,bin)
igg.bin=makeBin(igg,bin)

cons=host[which(host.bin%in%igg.bin),]
host.uniq=host[which(!host.bin%in%igg.bin),]
igg.uniq=igg[which(!igg.bin%in%host.bin),]

len=function(x){dim(x)[1]}

num.venn=draw.pairwise.venn(area1=len(host),area2=len(igg),cross.area=len(cons))

pdf(file.path(plotdir,"sv_venn.pdf"),width=3,height=3)
grid.draw(num.venn)
dev.off()

diff=do.call(rbind,list(host.uniq,igg.uniq))
diff.num=group_by(diff,type)%>%
    summarize(num=n())
diff.len=group_by(diff,type)%>%
    filter(!is.na(len),type!="TRA")%>%
    summarize(len=sum(as.numeric(len)))
igg.len=group_by(igg,type)%>%summarize(sum(as.numeric(len))

makeBars=function(df){
    ggplot(df,mapping=aes(x=x,y=y,color=x,fill=x))+
        geom_bar(stat="identity")+
        scale_y_log10()+
        theme_bw()+
        theme(legend.position="bottom")
}

df.len=rename(diff.len,x=type,y=len)
df.num=rename(diff.num,x=type,y=num)

g.len=makeBars(df.len) + labs(title="Total length of variants",x="Variant Type",y="Length")
g.num=makeBars(df.num) + labs(title="Total number of variants",x="Variant Type",y="Occurrences")

pdf(file.path(plotdir,"sv_bars.pdf"),width=6,height=4)
print(g.num)
print(g.len)
dev.off()


## ok let's try sliding window sv region plot
require(reshape2)

loc.cnames=c("chr","start","stop","type")
loc1=diff[,c(1:3,11)]%>%setNames(loc.cnames)
loc2=diff[,c(4:6,11)]%>%setNames(loc.cnames)
loc.bind=rbind(loc1,loc2)
loc=melt(loc.bind,id.vars=c("chr","type"),value.name="pos")[c("chr","pos","type")]
loc.tb=as.tibble(loc[order(loc$chr,loc$pos),])

loc.summary=loc.tb%>%
    filter(pos>0,type!="TRA") %>%
    group_by(chr) %>%
    summarize(num=n(),max=max(pos),min=min(pos),mean=mean(pos),range=max-min,freq=num/range)%>%
    arrange(desc(num))

chrom = loc.summary$chr[8]
## moving sum
bnum=10000
bp.chr=filter(loc.tb,chr==chrom) %>%
    mutate(bin=round(pos/bnum)*bnum) %>%
    group_by(bin,type) %>% summarize(count=n())
bp.chr.all=group_by(bp.chr,bin)%>%summarize(count=sum(count))

win=100000
wincount=tibble(start=seq(0,max(bp.chr$bin),1000),stop=start+win)
wincount$count=sapply(seq_along(wincount$start),function(x){
    sum(bp.chr.all$count[which(bp.chr.all$bin>wincount$start[x]&bp.chr.all$bin<=wincount$stop[x])])})

g1=ggplot(wincount,mapping=aes(x=start,y=count))+theme_bw()+
    geom_point(alpha=0.5,size=0.25)
g2=g1+geom_smooth(se=F,span=0.3,method="loess")
g3=ggplot(wincount,mapping=aes(x=start,y=count))+theme_bw()+
    geom_smooth(se=F,span=0.3,method="loess")


pdf(file.path(plotdir,"sv_win.pdf"),width=6,height=4,useDingbats=F)
print(g1)
print(g2)
print(g3)
dev.off()

##coverage
bamdir="/dilithium/Data/Nanopore/Analysis/171025_cho/bam/ngmlr"
covsuffix=".cov.txt"
hostfiles=system(paste0("ls ",bamdir,"/choNIHhost*",covsuffix),intern=T)
iggfiles=system(paste0("ls ",bamdir,"/choNIHIgG*",covsuffix),intern=T)
pooldir=file.path(bamdir,"../pooled")
hostpool=paste0(pooldir,"/choNIHhost",covsuffix)
iggpool=paste0(pooldir,"/choNIHIgG",covsuffix)

cnames=c("chr","pos","cov")
hostcov.list=lapply(seq_along(hostfiles),function(i){
    read_tsv(hostfiles[i],col_names=cnames)%>%
        mutate(replicate=rep(i))})
iggcov.list=lapply(seq_along(iggfiles),function(i){
    read_tsv(iggfiles[i],col_names=cnames)%>%
        mutate(replicate=rep(i))})
hostcov=do.call(rbind,hostcov.list)
iggcov=do.call(rbind,iggcov.list)
hostcov.pool=read_tsv(hostpool,col_names=cnames)%>%
    mutate(samp="host")
iggcov.pool=read_tsv(iggpool,col_names=cnames)%>%
    mutate(samp="igg")
cov.pool=rbind(hostcov.pool,iggcov.pool)

n=10000
hostcov.bin=mutate(hostcov,bin=round(pos/n)*n) %>%
    group_by(bin,replicate) %>% summarize(count=mean(cov))
iggcov.bin=mutate(iggcov,bin=round(pos/n)*n) %>%
    group_by(bin,replicate) %>% summarize(count=mean(cov))
hostcov.bin$replicate=factor(hostcov.bin$replicate)
iggcov.bin$replicate=factor(iggcov.bin$replicate)
cov.pool.bin=mutate(cov.pool,bin=round(pos/n)*n) %>%
    group_by(bin,samp) %>% summarize(count=mean(cov))

hostcov.bin$samp="host"
iggcov.bin$samp="igg"

## let's do moving average of the coverage
ylim=150
npoint=10000
cov.pool.point=mutate(cov.pool.bin,point=round(bin/npoint)*npoint)%>%
    group_by(point,samp) %>% summarize(num=mean(count)) %>%
    filter(num<=ylim)
host.smooth=filter(cov.pool.point,samp=="host")
igg.smooth=filter(cov.pool.point,samp=="igg")
require(zoo)
width=30
host.smooth$roll=rollapply(host.smooth$num,width=width,mean,fill=NA)
igg.smooth$roll=rollapply(igg.smooth$num,width=width,mean,fill=NA)
cov.smooth=rbind(host.smooth,igg.smooth)

cov.bin=rbind(hostcov.bin,iggcov.bin)

gcov=ggplot(cov.bin,mapping=aes(x=bin,y=count,color=samp,group=replicate))+theme_bw()+
    geom_point(alpha=0.3,size=0.1)+
    ylim(0,100)+
    geom_smooth(aes(group=samp,color=samp),se=F,span=0.3,method="loess") +
    labs(title="coverage along picr_49",x="Coordinate",y="Coverage") + theme(legend.position="bottom")
gcov.pool=ggplot(cov.pool.bin,mapping=aes(x=bin,y=count,color=samp,group=samp))+theme_bw()+
    geom_point(alpha=0.3,size=0.1)+
    ylim(0,ylim)+
    geom_line(data=cov.smooth,mapping=aes(x=point,y=roll,group=samp,color=samp),alpha=0.8)+
    #geom_smooth(,se=F,span=0.5,method="loess") +
    labs(title="coverage along picr_49",x="Coordinate",y="Coverage") + theme(legend.position="bottom")
gcov.host=ggplot(hostcov.bin,mapping=aes(x=bin,y=count,color=replicate,group=replicate))+theme_bw()+
    geom_point(alpha=0.5,size=0.25)+
    geom_smooth(se=F,span=0.3,method="loess") + labs(title="Host coverage",x="Coordinate",y="Coverage")+
    theme(legend.position="bottom")
gcov.igg=ggplot(iggcov.bin,mapping=aes(x=bin,y=count,color=replicate,group=replicate))+theme_bw()+
    geom_point(alpha=0.5,size=0.25)+
    ylim(0,ylim)+
    geom_smooth(se=F,span=0.3,method="loess") + labs(title="Producer coverage",x="Coordinate",y="Coverage")+
    theme(legend.position="bottom")

pdf(file.path(plotdir,"cov_win.pdf"),width=6,height=4,useDingbats=F)
print(gcov.pool)
#print(gcov.host)
#print(gcov.igg)
dev.off()
