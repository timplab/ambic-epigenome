source("/home/isac/Code/ilee/plot/ggplot_theme.R")
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)

root="/kyber/Data/NGS/projects/ambic/combined/data/atacseq"
outdir=file.path(root,"deseq2")
plotdir="~/Dropbox/Data/ambic/atac/plots"
counts.fp=file.path(root,"peaks/individual_peaks_counts.txt")

counts.tb = read_tsv(counts.fp,skip=1)
counts.tb = counts.tb[which(counts.tb$Chr!="MT"),] # remove MT
counts.tb = counts.tb[!duplicated(counts.tb[,-1]),] # remove duplicate peaks

## make pData
fpaths = colnames(counts.tb)[7:dim(counts.tb)[2]]
fnames = basename(fpaths)
samples = sapply(strsplit(fnames,"[.]"),"[[",1)
replicates = factor(as.numeric(str_sub(samples,-1)))
cells = factor(ifelse(grepl("CHOK1",samples),"CHOK1","CHOZN"))
types = factor(ifelse(grepl("ost",samples),"host","producer"))
media = factor(ifelse(grepl("Nogln",samples),"selective","unselective"))
stabilities = factor(ifelse(grepl("Unstable",samples),"Unstable","Stable"))
conditions = factor(sapply(strsplit(str_sub(samples,1,-2),"_"),"[[",1))
pd = tibble(sample=samples,cell=cells,type=types,
            media=media,stability=stabilities,
            condition=conditions,replicate=replicates)
chok1idx=grep("CHOK1",pd$sample)
choznidx=grep("CHOZN",pd$sample)

## convert to dds
cts = as.matrix(counts.tb[,7:dim(counts.tb)[2]],row.names="Geneid")
coords = counts.tb[,1:6]
coords.gr = GRanges(coords)
colnames(cts) = pd$sample
dds = DESeqDataSetFromMatrix(countData = cts,
                             colData = pd,
                             design = ~ condition,
                             rowRanges = coords.gr)

allgroups = c("cell","type","stability")

## prefilter based on count
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

## preprocess data
dds <- estimateSizeFactors(dds,locfunc = genefilter::shorth)
dds <- DESeq(dds,parallel=T)
colMeans(counts(dds,normalized=TRUE))

## save data
if (FALSE) {
    ddspath = file.path(outdir,"allsamples_dds.rds")
    saveRDS(dds,file=ddspath)
}

## QC
ntd <- normTransform(dds)
vsd <- vst(dds,blind=FALSE)

if (FALSE) {
    ## variance
    plotpath = file.path(plotdir,"atacseq_variance.pdf")
    pdf(plotpath,width=4,height=3,useDingbats=F)
    ## just chozn
    meanSdPlot(assay(vsd)) 
    ##meanSdPlot(assay(ntd))
    dev.off()
    ## using vsd for further analysis
}

if (TRUE) {
    ## pca comparison
    plotpath = file.path(plotdir,"atacseq_pca.pdf")
    pdf(plotpath,useDingbats=F)
    plotPCA(vsd,intgroup=allgroups)
    plotPCA(vsd,intgroup=c("cell"))
    plotPCA(vsd,intgroup=c("type"))
    plotPCA(vsd[,chok1idx],intgroup=allgroups)+
        labs(title="CHOK1")
    plotPCA(vsd[,choznidx],intgroup=allgroups)+
        labs(title="CHOZN")
    plotPCA(vsd[,grep("able",colData(dds)$sample)],intgroup=c("stability")) 
    plotPCA(vsd[,grep("able",colData(dds)$sample)],intgroup=c("media")) 
    dev.off()
}

if (TRUE) {
    ## dendrogram,heatmap?
    plotpath = file.path(plotdir,"atacseq_cluster.pdf")
    pdf(plotpath,useDingbats=F)
    datorder = order(rowMeans(counts(dds,normalized=TRUE)),
                     decreasing=TRUE)
    dists.chozn = dist(t(assay(vsd[,choznidx])))
    dists.chok1 = dist(t(assay(vsd[,chok1idx])))
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    for (dists in list(dists.chozn,dists.chok1)){
        distmat = as.matrix(dists)
        pd.sub = pd[pd$sample %in% rownames(distmat),]
        rownames(distmat) = paste(pd.sub$cell,pd.sub$type,pd.sub$stability,sep="-")
        colnames(distmat) = NULL
        pheatmap(distmat,clustering_distance_rows=dists,
                 clustering_distance_cols=dists,col=colors)
    }
    ann_df = colData(vsd)[,c("cell","type","stability")]
    pheatmap(assay(vsd)[datorder[1:10],],cluster_rows=FALSE,show_rownames=FALSE,
             cluster_cols=FALSE)# ,annotation_col=ann_df)
    dev.off()
}

## insertion?
insert.chozn = tibble(chrom="scaffold_16",start=15304000,end=start+1000)
insert.chok1 = tibble(chrom="scaffold_72",start=1634968,end=start+1000)
choznin.gr = GRanges(insert.chozn)
chok1in.gr = GRanges(insert.chok1)
coords.gr = GRanges(coords)
chozn.dist = as.tibble(distanceToNearest(coords.gr,choznin.gr))
chok1.dist = as.tibble(distanceToNearest(coords.gr,chok1in.gr))
chozn.idx = which.min(chozn.dist$distance)
chok1.idx = which.min(chok1.dist$distance)
indices = c(chozn.idx,chok1.idx)
plotpath = file.path(plotdir,"atacseq_insertflank_counts.pdf")
pdf(plotpath,height=3,width=5,useDingbats=F)
for (i in c(1,2)){
    idx = indices[i]
    counts = plotCounts(dds,gene=idx,
                        intgroup=c("cell","type","stability"),
                        returnData=T)

    if (i == 1){
        counts = counts[choznidx,]
    }else{
        counts = counts[chok1idx,]
    }
    counts$label = paste(counts$cell,counts$type,counts$stability,sep=":")
    g = ggplot(counts,aes(x=label,y=count,group=label,color=type))+
        geom_boxplot()+
        geom_jitter(height=0,width=0.1) +
        scale_y_log10()+
        coord_flip()+
        theme(legend.position="bottom")
    ##    theme(axis.text.x = element_text(angle=45,hjust=1,))
    print(g)
}
dev.off()

## comparison
res.chok1 = results(dds.chok1)
idxorder = order(res.chok1$padj)
res.chok1[idxorder,]
coords[idxorder,] %>% filter(Chr=="scaffold_93")
