library(tidyverse)
library(umap)
library(bsseq)

## Create a methylation umap from the BSseq and BSmooth R objects - updated for Kevin's use

plotdir <- "/home/kmcfarland/methylation_umap"
dir <- "/kyber/Data/Nanopore/projects/ambic/sigma/methylation/mfreq"
rawpath <- file.path(dir,"CHOZN.BSseq.rds")
bspath <- file.path(dir,"CHOZN.BSmooth.Rds")

# read data
#bsobj <- read.bismark(pd$fpath,colData = pd)
bsobj <- readRDS(bspath)
pd <- pData(bsobj) %>% as_tibble()
pd$rep <- as.factor(sapply(strsplit(pd$label,"_"),"[[",2)) #split by replicates
pd$day <- as.factor(sapply(strsplit(pd$sample,"Day"),"[[",2)) #Split by Days

# get methylation
meth <- getMeth(bsobj,type = "smooth")
meth <- as.matrix(meth)
meth <- na.omit(meth)
nrow(meth)
meth.t <- t(meth)

config <- umap.defaults
config$n_neighbors = 5
config$n_epochs = 500

## Get umap split only by cell line (+/- selection pressure comparison across days)
cell_lines <- list(c('CHOZNStableNogln','CHOZNStableGln'),c('CHOZNUnstableNogln','CHOZNUnstableGln'),c('CHOZNHost'))

fid <- c('Stable','Unstable','Host')
cnt <- 0
for (i in cell_lines){
cnt <- cnt+1
print(fid[cnt])
meth.umap <- umap(meth.t[pd$cell %in% i,],config)

umap.tb <- as.tibble(meth.umap$layout) %>% bind_cols(pd %>% filter(cell %in% i))

umap.out <- umap.tb %>% dplyr::select(-fpath)
outpath <- file.path(plotdir, paste('201119_methylation_',fid[cnt],'_umap.tsv',sep=''))
write_tsv(umap.out,outpath)
}


## Old file processing code
#fpaths <- system(paste("find",dir,"-name \"*txt.gz\""),intern=T)
#labs <- sapply(strsplit(basename(fpaths),"[.]"),"[[",1)
#samples <- sapply(strsplit(labs,"_"),"[[",1)
#cells <- sapply(strsplit(labs,"Day"),"[[",1)
#days <- as.factor(as.numeric(sapply(strsplit(samples,"Day"),"[[",2)))
#pd <- tibble(fpath = fpaths, label = labs, sample = samples, cell = cells, day = days)


# # By cell line and Gln selection

# for (i in unique(pd$cell)){
# print(i)
# meth.umap <- umap(meth.t[pd$cell %in% i,],config)
#
# umap.tb <- as.tibble(meth.umap$layout) %>% bind_cols(pd %>% filter(cell %in% i))
#
# umap.out <- umap.tb %>% dplyr::select(-fpath)
# outpath <- file.path(plotdir, paste('201119_methylation_',i,'_umap.tsv',sep=''))
# write_tsv(umap.out,outpath)




# # For cell line with combined GLN selection
# cell_lines <- list(c('CHOZNStableNogln','CHOZNStableGln'),c('CHOZNUnstableNogln','CHOZNUnstableGln'),c('CHOZNHost'))
#
# fid <- c('Stable','Unstable','Host')
# cnt <- 0
# for (i in cell_lines){
# print(i)
# fid[cnt]
# meth.umap <- umap(meth.t[pd$cell %in% i,],config)
#
# umap.tb <- as.tibble(meth.umap$layout) %>% bind_cols(pd %>% filter(cell %in% i))
#
# umap.out <- umap.tb %>% dplyr::select(-fpath)
# outpath <- file.path(plotdir, paste('201119_methylation_',fid[cnt],'_umap.tsv',sep=''))
# print(paste('Written to',outpath))
# write_tsv(umap.out,outpath)
# }
