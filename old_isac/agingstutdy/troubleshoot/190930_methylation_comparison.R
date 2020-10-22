#!/usr/bin/env Rscript
rm(list=ls());gc()
source("~/Code/ilee/plot/ilee_plot_utils.R")
library(bsseq)

dir <- "/kyber/Data/Nanopore/projects/ambic/sigma/methylation/subset_troubleshoot"
bedpaths <- system(paste("find",dir,"-name \"*bed.gz\""),intern=T)
runs <- sapply(strsplit(basename(bedpaths),"[.]"),"[[",1)
fpaths <- paste0(dir,"/",runs,".cpg.meth.freq.txt.gz")
labs <- sapply(strsplit(runs,"_"),"[[",1)
cells <- sapply(strsplit(labs,"Day"),"[[",1)
pd <- tibble(fpath = fpaths, run = runs, label = labs, cell = cells)
## choose one cell + one of each
#cell.select <- "CHOZNHost"
#cell.select <- "CHOZNStableGln"
#cell.select <- "CHOZNStableNogln"
#pd <- pd[
#  grepl(cell.select,cells) | 
#    runs %in% c("CHOZNHostDay30_1_PAD07455",
#      "CHOZNStableGlnDay30_2_PAD09863", 
#      "CHOZNStableNoglnDay30_3_PAD07512", 
#      "CHOZNUnstableGlnDay30_1_PAD05973",
#      "CHOZNUnstableNoglnDay60_1_PAD09584"), ] 

# read data
bsobj <- read.bismark(pd$fpath,colData = pd)

# get methylation
meth <- getMeth(bsobj,type = "raw")
meth <- na.omit(meth)
nrow(meth)
# excluding  weird sample
keepi <- !(grepl("FAH47313",pd$run)) 
pd <- pd[keepi,]
meth <- meth[,keepi]
ncol(meth)

# pca
meth.t <- t(meth)
pca.m <- prcomp(meth.t)
pca.tb <- as_tibble(pca.m$x[,1:6]) %>%
  bind_cols(pd)

# plot
g <- ggplot(pca.tb,aes(x = PC1, y = PC2, color = cell)) + geom_point()
plotter(g,"cho_pca")
g <- ggplot(pca.tb,aes(x = PC3, y = PC4, color = cell)) + geom_point()
plotter(g,"cho_pca2")

pca.tb[which(pca.tb$PC1 < 0 & pca.tb$cell == "CHOZNStableGln"),]$run
pca.tb[which(pca.tb$PC1 > 0 & pca.tb$cell == "CHOZNHost"),]$run

# outliers
pca.tb[which(pca.tb$PC2 > 40),]$run # FAH47313
pca.tb[which(pca.tb$PC3 < -40),]$run # FAH47656
pca.tb[which(pca.tb$PC4 < -40),]$run # FAH61777
