#!/usr/bin/Rscript
library(tidyverse)
library(reshape2)

#load
args=(commandArgs(TRUE))
csvfiles.string=args[1]
plotdir=args[2]
prefix=args[3]
print(args)
csvfiles=strsplit(csvfiles.string,"\n")[[1]]
##load files
cnames=c("rname", "rlength", "avgqual", "start.time", "end.time", "duration", "chnum", "mux", "bc")
ncols=count_fields(csvfiles[1],tokenizer_csv(),n_max=1)
stats.list=lapply(csvfiles,FUN=function(x){
    read_csv(x,col_names=cnames[1:ncols]) %>%
        mutate(rlength=as.numeric(rlength))})
##name list elements
sampnames=sapply(csvfiles,FUN=function(x){
n    strsplit(basename(x),"[.]")[[1]][1]})
names(stats.list)=sampnames
#cumSum?
stats.list=lapply(stats.list,
             FUN=function(x)
             {x %>%
                  arrange(rlength) %>%
                  mutate(cumlen=cumsum(rlength))})
n=1000 #binning number
sum.stats=lapply(stats.list,
           FUN=function(x)
           {x %>% summarize(i50=which.min(abs(cumlen-max(cumlen)/2)),
                            f50=i50/n(),
                            n50=rlength[i50],
                            b=round(n()/n), ##binning number
                            b50=i50/b/n, ## fraction of bases to n50
                            yield=max(cumlen),
                            nreads=n(),
                            med=median(rlength),
                            max=max(rlength),
                            len.lim=quantile(rlength,probs=.99))
           })

len.bin = mapply(
    function(x,y){
        stats.list[[x]] %>%
            group_by(bin=trunc(y$b:(n()+y$b-1)/y$b))  %>%
                summarize(cumlen=mean(cumlen),rlen=mean(rlength))  %>%
                transmute(samp=x,fread=bin/n,cumf=cumlen/max(stats.list[[x]]$cumlen),rlen=rlen)
    },
    x=names(stats.list),y=sum.stats,SIMPLIFY=F) %>%
    bind_rows

##plot
len.df=lapply(
    names(stats.list),
    function(x){
        stats=stats.list[[x]]
        data.frame(samp=x,rlength=stats$rlength)
    }) %>%
    bind_rows
##specific for cho
rep=sapply(strsplit(len.df$samp,"_choNIHhost"),function(x) x[2])
len.df$rep=rep
###
sum.df=sum.stats %>%
    bind_rows %>%
    mutate(samp=names(sum.stats))

n50=melt(sapply(sum.stats,function(x) x$n50)) %>%
    mutate(samp=names(sum.stats))
len.lim=max(sum.df$len.lim)
p.ecdf=ggplot(len.bin,aes(x=rlen,group=samp,color=samp))+theme_bw() +
    stat_ecdf()+
    geom_vline(data=n50,aes(xintercept=value,color=samp),
               linetype="dashed",size=0.5,alpha=0.5) +
    xlim(0,len.lim)+
    labs(title="Read length Empirical Cumulative Distribution", y="Fraction of Reads",x="Read Length")
p.cum=ggplot(len.bin,aes(x=rlen,y=cumf,group=samp,color=samp))+theme_bw() +
    geom_line()+
    geom_vline(data=sum.df,aes(xintercept=n50,color=samp),
               linetype="dashed",size=0.5,alpha=0.5) +
    xlim(0,len.lim)+
    labs(title="Cumulative Plot of Read Length", y="Fraction of Total Yield",x="ReadsLength")


fname=file.path(plotdir,paste0(prefix,"_plots.pdf"))
pdf(file=fname,width=9,height=6)
print(p.ecdf)
print(p.cum)
dev.off()
