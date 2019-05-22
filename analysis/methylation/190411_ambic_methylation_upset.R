source("/home/isac/Code/ilee/plot/ggplot_theme.R")
library(GenomicRanges)
library(UpSetR)
library(ComplexHeatmap)
library(ensembldb)

root="/kyber/Data/Nanopore/projects/ambic/combined/rds"
outdir=root
plotdir="~/Dropbox/Data/ambic/methylation/plots"
gffpath = "/mithril/Data/NGS/Reference/cho/chok1/annotation/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.96.gff3"

allgroups = c("cell","type","stability")
## load data
if (TRUE) {
    setwd(root)
    stablegln = readRDS("DMR.Stable.Gln.rds")
    stablenogln = readRDS("DMR.Stable.Nogln.rds")
    unstablegln = readRDS("DMR.Unstable.Gln.rds")
    chok1 = readRDS("DMR.K1.DMRSeq.rds")
    dmr.list = list(chok1,stablegln,stablenogln,unstablegln)
    names(dmr.list) = c("CHOK1IgG","CHOZNStableGlnDay0",
                        "CHOZNStableNoglnDay0","CHOZNUnstableGlnDay0")
}
## subsetting by qvalue criterion
if (TRUE){
    alpha = 0.2
    sig.list = lapply(dmr.list,function(x){
      x[which(x$qval<=alpha)]
    })
}

## genes
gff = rtracklayer::import.gff(gffpath)
gff.genes = gff[grep("gene",gff$ID)]
genes.prom = promoters(gff.genes)
## enrichment in promoters?
promwidth = sum(width(genes.prom))/1e9
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
  
plotpath = file.path(plotdir,"methylation_dmr_promoter_enrichment.pdf")
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
plotpath = file.path(plotdir,"dmr_upset.pdf")
pdf(plotpath,height=3,width=6,useDingbats=F)
#upset(fromList(sig.ind),order.by="freq")
UpSet(sig.mat,comb_order = order(comb_size(sig.mat),decreasing=T))
dev.off()


