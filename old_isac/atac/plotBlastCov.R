#!/usr/bin/Rscript
library(GenomicRanges)
library(ggplot2)
# define the inputs here
samp="choATACIgG1"
indir="/atium/Data/NGS/Aligned/170530_choATAC/blast"
intxt=file.path(indir,paste0(samp,"_blastn.txt"))
mitoname="gi|91176202|ref|NC_007936.1|"
outdir="/home/isac/Dropbox/Data/ambic/atac/plot"
outpath=file.path(outdir,paste0(samp,"_mito.pdf"))

# read the blastn text output
blast=read.table(intxt,stringsAsFactors=F)
blast.df=data.frame(blast[,c(2,9,10)])
colnames(blast.df)=c("seqnames","start","end")

mito.df=blast.df[which(blast.df$seqnames==mitoname),]
mito.order=order(mito.df$start)
mito.df=mito.df[mito.order,]

#flip reverse alignments
widths=mito.df$end-mito.df$start
mito.df[which(widths<0),c(2,3)]=mito.df[which(widths<0),c(3,2)]

# make granges object
mito.gr=GRanges(seqnames=mito.df$seqnames,ranges=IRanges(start=mito.df$start,end=mito.df$end))
mito.cov=as.vector(coverage(mito.gr)[[1]])

cov.df=data.frame(pos=seq(1,length(mito.cov)),cov=mito.cov)

g = ggplot(cov.df,aes(x=pos,y=cov))+geom_point()

pdf(outpath)
print(g)
dev.off()
