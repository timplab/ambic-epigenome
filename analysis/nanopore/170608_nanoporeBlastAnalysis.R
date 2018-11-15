library(GenomicRanges)
library(GenomicAlignments)
library(ggplot2)
library(reshape2)

plotdir="~/Dropbox/Data/ambic/nanopore/plots"
outdir="/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/analysis"
bamdir="/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/oldgenome/"
bampath="/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/oldgenome/170608_choIgGNIH_hits.sorted.bam"
bampath.gmap="/atium/Data/Nanopore/Analysis/170608_cho/rawhit/hits.graphmap.sorted.bam"
adir="/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/analysis"

##bam inputs
param=ScanBamParam(what="mapq")
bambp.names=read.table(file.path(adir,"bpNames.txt"),stringsAsFactors=F)[,1]
bam.all=readGAlignments(bampath,param=param)
bam.gmap.all=readGAlignments(bampath.gmap,param=param)
bam=bam.all[which(mcols(bam.all)$mapq>20)]
bam.gmap=bam.gmap.all[which(mcols(bam.gmap.all)$mapq>20)]
bam.names.all=readGAlignments(bampath,use.names=T,param=param)
bam.names=bam.names.all[which(mcols(bam.all)$mapq>20)]
##bam quality plots

bwa.q.df = data.frame(mapq=mcols(bam.all)$mapq)
gmap.q.df = data.frame(mapq=mcols(bam.gmap.all)$mapq)
g.bwaq = ggplot(bwa.q.df,aes(x=mapq))+geom_density()+
    theme_bw() + labs(x="map quality",y="density",title="bwa alignment map quality")
g.gmapq = ggplot(gmap.q.df,aes(x=mapq))+geom_density()+
    theme_bw() + labs(x="map quality",y="density",title="bwa alignment map quality")

pdf(file.path(plotdir,"mapQuality.pdf"))
print(g.bwaq)
print(g.gmapq)
dev.off()

##granges for the breakpoints
bp.gr = GRanges(seqnames=c("KE379390","KE379503","KE379503","KE382060"),ranges=IRanges(start=c(131147,485274,6631618,5052404),width=1))
##subsetting just one of the regions
bpone.gr=bp.gr[4]
bpone.up=resize(bpone.gr,width=1000,fix="end")
end(bpone.up)=end(bpone.up)-100
bpone.dn=resize(bpone.gr,width=100000,fix="start")
start(bpone.dn)=start(bpone.dn)+100
chrname=as.character(seqnames(bpone.gr))
pname="4.1DHFR_VRC01WTG1M3_DGV"
##find overlap bams

upidx=queryHits(findOverlaps(bam.names,bpone.up))
dnidx=queryHits(findOverlaps(bam.names,bpone.dn))

bpidx=queryHits(findOverlaps(bam.names,bpone.gr))
names.ovl=names(bam.names[bpidx])
allbpidx=which(names(bam.names) %in% names.ovl)
bam.bp=bam.names[bpidx]
bam.ins=bam.names[allbpidx]

names.max=names(bam.ins[which(mcols(bam.ins)$mapq==60)])
max.freq=as.data.frame(table(names.max))
max.freq=max.freq[order(max.freq$Freq,decreasing=T),]
max.name=as.character(max.freq$names.max[1])
maxbam=bam.names[which(names(bam.names)==max.name)]
sum(width(maxbam))
max.names=as.character(max.freq$names.max[1:3])
write.table(x=max.names,file=file.path(outdir,"egbams.txt"),quote=F,row.names=F,col.names=F)

bam.df=data.frame(chr=character(),start=numeric(),end=numeric(),ymin=numeric(),ymax=numeric())
for (i in seq(length(names.ovl))){
    n = names.ovl[i]
    b = bam.ins[which(names(bam.ins)==n)]
    df = data.frame(chr=seqnames(b),start=start(b),end=end(b),ymin=i-1,ymax=i)
    plen = sum(width(b[which(seqnames(b)==pname),]))
    plen.df = data.frame(chr="plasmid",start=0,end=plen,ymin=i-1,ymax=i)
    bam.df=rbind(bam.df,df,plen.df)
}


##plot the overlap bams using rect
#p1 is plot for genome part
p1.df=bam.df[which(bam.df$chr==chrname),]
p1.df=p1.df[which(p1.df$start<end(bpone.gr)),]
#p2 is plot for plasmid part
p2.df=bam.df[which(bam.df$chr==pname),]
#p3 is the plot of plasmid length
p3.df=bam.df[which(bam.df$chr=="plasmid"),]
g1 = ggplot(p1.df,aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax))+geom_rect()+
    theme_bw() + labs(x="position",y="read",title="pileup on cho genome")
g2 = ggplot(p2.df,aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax))+geom_rect()+
    theme_bw() + labs(x="position",y="read",title="pileup on plasmid")
g3 = ggplot(p3.df,aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax))+geom_rect()+
    theme_bw() + labs(x="length",y="read",title="length of plasmid alignments")
pdf(file.path(plotdir,"insertionBP_bwa.pdf"))
print(g1)
print(g2)
print(g3)
dev.off()

## plasmid related plots
bam.bwa.p=bam[which(seqnames(bam)==pname)]
bwa.pcov=coverage(bam.bwa.p)[[pname]]
bwa.pcov.num=as.numeric(bwa.pcov)
bwa.pcov.df=data.frame(start=seq(length(bwa.pcov.num)),cov=bwa.pcov.num)
bwa.df=as.data.frame(bam.bwa.p)[c("seqnames","start","end","width")]
bwa.df=cbind(bwa.df,data.frame(ymin=seq(0,length(bam.bwa.p)-1),ymax=seq(length(bam.bwa.p))))
g.bwa.cov=ggplot(bwa.pcov.df,aes(x=start,y=cov)) + geom_point(size=0.5) +
    theme_bw()+ labs(x="position",y="coverage",title="plasmid alignment coverage")
g.bwa.a=ggplot(bwa.df,aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax))+geom_rect()+
    theme_bw()+ labs(x="position",y="read",title="plasmid alignment pileup")
g.bwa.den = ggplot(bwa.df,aes(x=width))+geom_density()+geom_vline(xintercept=median(wp),color="blue")+
    theme_bw() + labs(x="alignment length",y="density",title="plasmid alignment length density plot")
pdf(file.path(plotdir,"plasmid_alignment_bwa.pdf"))
print(g.bwa.cov)
print(g.bwa.a)
print(g.bwa.den)
dev.off()

## graphmap specific plots
bam.p=bam.gmap[which(seqnames(bam.gmap)==pname)]
p.cov=coverage(bam.p)[[pname]]
p.cov.num=as.numeric(p.cov)
p.cov.df=data.frame(start=seq(length(p.cov.num)),cov=p.cov.num)
p.df=as.data.frame(bam.p)
p.df=p.df[c("seqnames","start","end","width")]
p.df=cbind(p.df,data.frame(ymin=seq(0,length(bam.p)-1),ymax=seq(length(bam.p))))
g.cov=ggplot(p.cov.df,aes(x=start,y=cov)) + geom_point(size=0.5) +
    theme_bw() + labs(x="position",y="coverage",title="plasmid alignment coverage")
g.a=ggplot(p.df,aes(xmin=start,xmax=end,ymin=ymin,ymax=ymax))+geom_rect()+
    theme_bw() + labs(x="position",y="read",title="plasmid alignment pileup")
g.den = ggplot(p.df,aes(x=width))+geom_density()+geom_vline(xintercept=median(wp),color="blue")+
    theme_bw() + labs(x="alignment length",y="density",title="plasmid alignment length density plot")
pdf(file.path(plotdir,"plasmid_alignment_graphmap.pdf"))
print(g.cov)
print(g.a)
print(g.den)
dev.off()
