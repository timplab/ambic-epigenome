source("/home/isac/Code/ilee/plot/ggplot_theme.R")
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)

root="/kyber/Data/NGS/projects/ambic/combined/data/atacseq"
outdir=file.path(root,"deseq2")
plotdir="~/Dropbox/Data/ambic/atac/plots"

allgroups = c("cell","type","stability")
## load data
if (TRUE) {
    ddspath = file.path(outdir,"allsamples_dds.rds")
    dds = readRDS(file=ddspath)
}
pd = colData(dds)
chok1idx=grep("CHOK1",pd$sample)
choznidx=grep("CHOZN",pd$sample)
chozn.conditions = unique(as.character(pd$condition[grep("ln",as.character(pd$condition))]))
chok1.conditions = unique(as.character(pd$condition[grep("IgG",as.character(pd$condition))]))
chozn.host = unique(as.character(pd$condition[grep("CHOZNHost",as.character(pd$condition))]))
chok1.host = unique(as.character(pd$condition[grep("CHOK1host",as.character(pd$condition))]))

## separate chok1 vs chozn
dds.chok1 = dds[,chok1idx]
dds.chozn = dds[,choznidx]

##  drop levels and designate reference
dds.chozn$condition <- relevel(droplevels(dds.chozn$condition),
                               ref = chozn.host)
dds.chok1$condition <-  relevel(droplevels(dds.chok1$condition),
                                ref = chok1.host)

## redo prefilter and normalization
dds.list = list(dds.chok1,dds.chozn)
dds.list = lapply(dds.list,function(x){
    keep = rowSums(counts(x)) >= 5
    x[keep,]
})

## comparisons
dds.list = lapply(dds.list,function(x){
    DESeq(x,parallel=T)
})

## overlaps?
dds.chok1 = dds.list[[1]]
dds.chozn = dds.list[[2]]

alpha = 0.1
sig.list = lapply(chozn.conditions,function(x){
    contrast=c("condition",x,chozn.host)
    res = results(dds.chozn,contrast)
    rowRanges(dds.chozn)[which(res$padj<alpha)]
})
names(sig.list) = chozn.conditions
chok1.res = results(dds.chok1)
chok1.sig = rowRanges(dds.chok1)[which(chok1.res$padj<alpha)]
sig.list[[chok1.conditions]] = chok1.sig

coords.list = lapply(seq_along(sig.list),function(i){
    name = names(sig.list)[i]
    x = as.tibble(sig.list[[name]])  %>%
        select(Geneid,seqnames,start,end) %>%
        mutate(sample=name)
})
sig.coords = do.call(rbind,coords.list)
sig.spread = sig.coords %>%
    spread(sample,Geneid)

## suppvec
dat.supp = sig.spread[,4:dim(sig.spread)[2]]
samples = colnames(dat.supp)
suppvec = apply(ifelse(is.na(dat.supp),0,1),
                1,paste,collapse="")
suppnum = group_by(as.tibble(suppvec),value)%>%
    summarize(n = n()) %>%
    arrange(desc(n))




