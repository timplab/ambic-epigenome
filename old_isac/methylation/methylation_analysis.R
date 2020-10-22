library(tidyverse)
library(ggdendro)
library(seriation)


workdir="/atium/Data/Nanopore/Analysis/170608_cho/rawhit/"
plotdir="~/Dropbox/Data/ambic/170608_hitreads"


freq=read_tsv(file.path(workdir, "frequencies_ngmlr.tsv"))

plasmid.name="4.1DHFR_VRC01WTG1M3_DGV"

plasmid.freq=freq %>%
    filter(chromosome==plasmid.name)

pdf(file.path(plotdir, "try.pdf"))

ggplot(plasmid.freq, aes(x=start, y=methylated_frequency))+geom_line()+theme_classic()

dev.off()

phasedat=read_tsv(file.path(workdir, "methylation_bwa.phase.tsv")) %>%
    filter(chr==plasmid.name) %>%
    mutate(tempname=unlist(lapply(strsplit(readix, split=":"), function(x) {x[1]})))




##ok - what to do.  Need to take all CGs that exist, make an array of 0s, reshape, etc.

phasedat$masked=phasedat$meth
phasedat$masked[phasedat$masked==0]=NA
phasedat$masked[phasedat$masked==-1]=0

sv.reads=read_tsv(file.path(workdir, "bpinfo.txt"), col_names=F)
colnames(sv.reads)=c("Insertion", "template", "file")
sv.reads$template=unlist(strsplit(sv.reads$template,split=":"))[rep(c(T,F,F))]

highmapq=read_tsv(file.path(plotdir, "egbams.txt"), col_names=F)


phasedat.sv=phasedat %>%
    filter(tempname %in% sv.reads$template) %>%
    mutate(insert=sv.reads$Insertion[pmatch(tempname, sv.reads$template, duplicates.ok=TRUE)])



phasedat.high=phasedat %>%
    filter(tempname %in% highmapq$X1)

insert.meth=phasedat.sv %>%
    group_by(insert) %>%
    summarize(avg.meth=mean(masked, na.rm=T), numcalled=sum(!is.na(masked)), num.reads=length(unique(readix)))

insert.meth=ddply(phasedat.sv,.(insert),function(x) c(avg.meth=mean(x$masked,na.rm=T),numcalled=sum(!is.na(x$masked)),num.reads=length(unique(x$readix))))

group.meth=na.omit(ddply(phasedat.sv,.(insert,start),function(x) mean(x$masked,na.rm=T,)))

g = ggplot(data=group.meth,mapping=aes(x=start,y=V1,group=insert,color=insert))
gfull = g + geom_point(alpha=0.5)+
    labs(x="Plasmid Coordinate",y="Methylation")+
    theme_bw()
gzoom = g + geom_point(alpha=0.5,size=0.7)+geom_smooth(se=F) +
    scale_x_continuous(limits=c(11000,12784))+
    labs(x="Plasmid Coordinate",y="Methylation")+
    theme_bw()

pdf(file.path(plotdir, "insertOverlay.pdf"),width=8,height=5)
print(gfull)
print(gzoom)
dev.off()
