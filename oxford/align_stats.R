#!/usr/bin/Rscript
# arguments are
# 1) full path to the summary file
# 2) full path (including .txt ext) to the output file
library(ggplot2)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)
fpath=args[1]
outpath=args[2]
print(fpath)
fname=basename(fpath)
base=strsplit(fname,"[.]")[[1]][1]

info=read.table(gzfile(fpath),sep=",",header=F,stringsAsFactors=F)
lengths=as.numeric(info[order(info[,2]),2])
print(max(lengths))
qual=info[,3]
stime=info[,4]
etime=info[,5]

qnum=4
qreads=length(lengths)/qnum
qseq=floor(seq(qreads,length(lengths),by=qreads))
qseq[length(qseq)]=length(lengths)
quantiles=lengths[qseq]
lenseq=c(2000,5000,10000,50000,100000)
numlen=numeric()
for (x in lenseq){
    numlen=c(numlen,(length(lengths)-which(lengths>x)[1]))
}
numwithin=which(lengths>lenseq[1])[1]
for (i in seq(length(numlen)-1)){
    numwithin=c(numwithin,numlen[i]-numlen[i+1])
}
numwithin = c(numwithin,numlen[length(numlen)])
num.df = t(data.frame(numwithin))
colnames(num.df)=c(0,lenseq)
quant.df = t(data.frame(quantiles))
colnames(quant.df)=seq(qnum)

stat.df=data.frame(baseGb=sum(lengths)/1000000000,numread=length(lengths),med=median(lengths),mean=mean(lengths),quant.df,num.df)
colnames(stat.df)=c("baseGb","numread","med","mean",seq(qnum),0,lenseq)

print(paste0("writing summary to ",outpath))
write.table(x=stat.df,file=outpath,quote=F,sep="\t",row.names=F)
