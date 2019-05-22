source("/home/isac/Code/ilee/plot/ggplot_theme.R")
library(GenomicRanges)
library(DESeq2)
library(ensembldb)
library(UpSetR)
library(ComplexHeatmap)

root="/kyber/Data/NGS/projects/ambic/combined/data/atacseq"
outdir=file.path(root,"deseq2")
plotdir="~/Dropbox/Data/ambic/atac/plots"
gffpath = "/mithril/Data/NGS/Reference/cho/chok1/annotation/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.96.gff3"

allgroups = c("cell","type","stability")
## load data
if (FALSE) {
    ddspath = file.path(outdir,"allsamples_dds.rds")
    dds = readRDS(file=ddspath)
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
    outpath = file.path(outdir,"dds_comparisons.rds")
    saveRDS(dds.list,file=outpath)
}
  
if (FALSE){
    ddspath = file.path(outdir,"dds_comparisons.rds")
    dds.list = readRDS(file=ddspath)
    dds.chok1 = dds.list[[1]]
    dds.chozn = dds.list[[2]]
    
    alpha = 0.05
    sig.list = lapply(chozn.conditions,function(x){
        contrast=c("condition",x,chozn.host)
        res = results(dds.chozn,contrast)
        rowRanges(dds.chozn)[which(res$padj<alpha)]
    })
    names(sig.list) = chozn.conditions
    chok1.res = results(dds.chok1)
    chok1.sig = rowRanges(dds.chok1)[which(chok1.res$padj<alpha)]
    sig.list[[chok1.conditions]] = chok1.sig
    sig.list = sig.list[c(5,1,2,3,4)]
    outpath = file.path(outdir,"comparison_sigdiff_ranges.rds")
    saveRDS(sig.list,file=outpath)
}

if (TRUE){
    sigpath = file.path(outdir,"comparison_sigdiff_ranges.rds")
    sig.list = readRDS(file=sigpath)
}

## genes
gff = import.gff(gffpath)
gff.genes = gff[grep("gene",gff$ID)]
w = 1000
genes.prom = promoters(gff.genes)
genes.ranges = resize(gff.genes,width = width(gff.genes)+w)
## enrichment in promoters?
promwidth = sum(width(genes.prom))/1000000000
genomesize = 2.4
sigprom.list = lapply(sig.list,
                      subsetByOverlaps,genes.prom)
promfrac = sapply(seq_along(sig.list),function(i){
  length(sigprom.list[[i]])/length(sig.list[[i]])
})
promrich = promfrac/(promwidth/genomesize)
promrich.tb = tibble(sample = names(sig.list),
                     enrichment = promrich)
## color by selection
selection = rep("No",dim(promrich.tb)[1])
selection[grep("Nogln",names(sig.list))] = "Yes"
promrich.tb$Selection = selection
## bar chart
g = ggplot(promrich.tb,aes(x=sample,y=enrichment,fill=Selection)) +
  geom_histogram(stat="identity") +
  theme(axis.text.x = element_text(angle=45,vjust=.5))
  
plotpath = file.path(plotdir,"atacseq_diffacc_promoter_enrichment.pdf")
pdf(plotpath,height=4,width=4,useDingbats=F)
print(g)
dev.off()

## upset plot
## make overlaps based on feature, not nucleotide base
gr.all = GRanges()
for (gr in sig.list){
  gr.all = c(gr.all,gr)
}
gr.all = unique(gr.all)
sig.ind = lapply(sig.list,function(x){
  subjectHits(findOverlaps(x,gr.all))
})

## make the combination matrix
sig.mat = make_comb_mat(sig.ind)
comb_size(sig.mat)[order(comb_size(sig.mat),decreasing=T)]
set_size(sig.mat)

## plot
plotpath = file.path(plotdir,"atacseq_upset.pdf")
pdf(plotpath,height=3,width=6,useDingbats=F)
#upset(fromList(sig.ind),order.by="freq")
UpSet(sig.mat,comb_order = order(comb_size(sig.mat),decreasing=T))
dev.off()


