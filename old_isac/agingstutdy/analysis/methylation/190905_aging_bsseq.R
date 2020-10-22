rm(list=ls());gc()
source("/home/isac/Code/ilee/plot/ggplot_theme.R")
library(bsseq)

# data info ----
root <- "/home/isac/Data/ambic/pooled_rep/mfreq"
regpath <- "/mithril/Data/NGS/Reference/cho/picr_ensembl/annotation/cho_picr_genes.bed"
days <- factor(c(0,30,60,90))
reps <- factor(c(1,2,3))
samples <- c("Host","StableGln","StableNogln","UnstableGln","UnstableNogln")

pd <- tibble(sample = rep(samples,length(days)*length(reps)),
             day = rep(rep(days,each = length(samples)),each = length(reps)),
             replicate = rep(rep(reps,each = length(samples)),length(days)),
             name = paste0("CHOZN",sample,"Day",day,"_",replicate),
             fp = paste0(root,"/",name,".cpg.meth.freq.txt.gz"))
bspath <- file.path(root,"choSigmaAgingStudy_BSseq.Rds")

# read as bsseq object ----
if ( ! file.exists(bspath)){
  bsobj <- read.bismark(pd$fp,colData = pd)
  saveRDS(bsobj,bspath)
} else {
  bsobj <- readRDS(bspath)
}

# read gene annotations ----
regs <- read_tsv(regpath,col_names = c("chrom","start","end","name","score","strand","id")) %>%
  mutate(start = start + 1)
regs.gr <- GRanges(regs)
proms <- promoters(regs.gr)

# data in regions ----
meth.regs <- getMeth(bsobj,proms,type = "raw",what = "perRegion")
meth.regs <- as_tibble(meth.regs)
names(meth.regs) <- pData(bsobj)$name
meth.regs$idx <- seq_along(regs.gr)

# gather ----
meth.gather <- meth.regs %>%
  gather(name,freq,-idx) %>%
  mutate(sample = rep(pd$sample,each = length(regs.gr)),
         day = rep(pd$day,each = length(regs.gr)),
         replicate = rep(pd$replicate,each = length(regs.gr)))

# collapse replicates ----
meth.avg <- meth.gather %>%
  group_by(sample,day,idx) %>%
  summarize(freq = mean(freq)) %>%
  na.omit()

# first density plots? ----
ggplot(meth.avg,aes(x = freq, group = day, color = day)) +
  facet_wrap(~sample) +
  geom_density()

ggplot(meth.avg,aes(y = freq, color = day, x = sample)) + 
  geom_violin(aes(fill = day),alpha = 0.3,width = 0.7) + geom_boxplot(width = 0.7,alpha = 0) +
  labs( x = "Sample", y = "Methylation", title = "Distributions of promoter methylation")
  
# insertion is near :
ins <- GRanges(tibble(chrom="RAZU01001824.1",start = 199405, end = 200000))
nearest_gene <- proms[nearest(ins,proms)]
