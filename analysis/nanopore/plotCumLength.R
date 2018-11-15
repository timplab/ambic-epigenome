#!/usr/bin/Rscript
library(tidyverse)

#load
args=(commandArgs(TRUE))
csvfile=args[1]
plotdir=args[2]
prefix=args[3]
cnames=c("rname", "rlength", "avgqual", "start.time", "end.time", "duration", "chnum", "mux", "bc")
ncols=count_fields(csvfile,tokenizer_csv(),n_max=1)
stats=read_csv(csvfile,col_names=cnames[1:ncols]) %>%
  mutate(rlength=as.numeric(rlength))
#cumSum?
stats=stats %>%
    arrange(stats$rlength) %>%
    mutate(cumlen=cumsum(rlength))
n50=as.numeric(stats %>% summarize(rlength[which.min(abs(cumlen-max(cumlen)/2))]))
n50ind=which(stats$rlength==n50)[1]
overn50num=length(stats$rlength)-n50ind
##binning
n=1000
b = round(length(stats$cumlen)/n)
bn50ind=round(n50ind/b)
f50=bn50ind/n
binlen = stats %>%
    group_by(x=trunc(b:(n()+b-1)/b)) %>%
    summarize(mean=mean(cumlen)) %>%
    transmute(f=x/n,flen=mean/max(stats$cumlen))
#plot
gp=ggplot(binlen,aes(x=f,y=flen))+theme_bw() +
    geom_line()+geom_vline(xintercept=f50,linetype="dashed",show.legend=TRUE)+
    labs(y="Fraction of total length",
         x="Fraction of reads")
fname=file.path(plotdir,paste0(prefix,"_cumlen.pdf"))
pdf(file=fname)
print(gp)
dev.off()
