#!/usr/bin/Rscript
library(tidyverse)
library(bsseq)
library(getopt)
dir=dirname(get_Rscript_filename())
if (is.na(NA)) dir = "."
source(file.path(dir,"../../../methylation/methylation_plot_utils.R"))
cores = ceiling(detectCores()/2)

pdpath = "/dilithium/Data/chosigma_meth_samples.txt"
outpath = "/dilithium/Data/tmp.rds"
mfreqdir = "/dilithium/Data/mfreq"
regspath = "/dilithium/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.fa.bins.shuffle20k.bed"
regspath = "/dilithium/Data/NGS/Reference/cho/chok1/genome/sigmaIgG.bed"
plotdir = "/dilithium/Data/plots"

pd = read_tsv(pdpath)

pd$filepath = file.path(mfreqdir,paste0(pd$samples,".chok1gshd_sigmaIgG.methfreq.txt.gz"))

# get frequencies
mfreq.list = mclapply(pd$filepath,mc.cores=cores,function(x){
    tabix_mfreq(x,regspath,cov=0,trinuc_exclude=NULL)
})
# add label column
for (i in seq_along(pd$samples)){
    mfreq.list[[i]]$samp = pd$samples[i]
}

# turn into bsseq object
mfreq.all = do.call(rbind,mfreq.list)
pd$samp = pd$samples
mfreq.bsseq = mfreqToBsseq(mfreq.all)
pData(mfreq.bsseq) = pd

# smooth
bs.cov = getCoverage(mfreq.bsseq,type="Cov",what="perBase")
keepi = which(rowSums(bs.cov>1)>dim(mfreq.bsseq)[2]/4*3)
bs.keep = mfreq.bsseq[keepi,]
bpparam = MulticoreParam(workers = cores,progressbar = TRUE)
bs.fit = BSmooth(bs.keep,BPPARAM=bpparam,verbose=TRUE,ns=20,h=500)

# PCA?
meth.smooth = getMeth(bs.fit,type="smooth",what="perBase")
meth = na.omit(meth.smooth)
meth.t = t(meth)
meth.t = meth.t[which(pd$samples!="choSigmaUnstableD0Glutrep3" | pd$samples != "choSigmaUnstableD0NoGlutrep3" ),]

meth.pca = prcomp(meth.t,center=TRUE,retx=TRUE,scale.=TRUE)

meth.x = bind_cols(pd,data.frame(meth.pca$x[,1:6]))
meth.x = bind_cols(pd[which(pd$samples!="choSigmaUnstableD0Glutrep3" | pd$samples != "choSigmaUnstableD0NoGlutrep3"),],data.frame(meth.pca$x[,1:6]))


# plot
g = ggplot(meth.x,aes(x=PC1,y=PC2,colour=cell))+geom_point()+theme_bw()
g2 = ggplot(meth.x,aes(x=PC1,y=PC2,colour=factor(day)))+geom_point()+theme_bw()
g3 = ggplot(meth.x,aes(x=PC1,y=PC2,colour=media))+geom_point()+theme_bw()
g4 = ggplot(meth.x,aes(x=PC1,y=PC2,colour=factor(replicate)))+geom_point()+theme_bw()

pdf(file.path(plotdir,"methylationPcaBiplot.pdf"),useDingbats=F)
print(g)
print(g2)
print(g3)
print(g4)
dev.off()
