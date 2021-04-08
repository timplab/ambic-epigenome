#!Rscript

## If a package is installed, it will be loaded. If any
## are not, the missing package(s) will be installed
## from CRAN and then loaded.

## First specify the packages of interest
# packages = c("tidyverse", "optparse",
#              "cluster", "DNAcopy","parallel")
#
# ## Now load or install&load all
# package.check <- lapply(
#   packages,
#   FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#       install.packages(x, repos = "http://cran.us.r-project.org", dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#   }
# )
# plot average curves


## FUNCTION : THIS SCRIPT TAKES PREVIOUSLY MADE COUNT FILES (binBedcounts.py) AND RUNS CNV ANALYSIS IN COMPARISON TO A GIVEN GENOME VIA .FAIDX FILE FOR A GIVEN WINDOW / BINWIDTH
library(optparse)
library(tidyverse)
library(cluster)
library(DNAcopy)
library(parallel)

parseArgs = function(argsin){
    arglist=list(
        make_option(c("-o","--output"),type="character",default=NULL,
                    help="output pdf file path",metavar="/path/to/out"),
        make_option(c("-w","--window"),type="numeric",default=50000,
                    help="binning window - should match binned cov window"),
        make_option(c("-f","--faidx"),type="character",default=NULL,
                    help="faidx for genomewide plot",metavar="/path/to/faidx")
    )
    argparser=OptionParser(option_list=arglist)
    args=parse_args(argparser,args=argsin,positional_arguments=c(1,100))
    ##Check if there are any variables, and if not, show help
    if (is.null(args$options$output)){
        print_help(argparser)
        print(args)
        stop("output must be provided",call.=FALSE)
    }
    args
}

readcounts <- function(fpath){
    cnames=c("chrom","start","end","count","perbp")
    dat=read_tsv(fpath,col_names=cnames)
    dat
}
trimCov <- function(dat,faidx,binwin){
    keepchrom=faidx$chrom[1:40]
    dat.filt=dat[which(dat$chrom %in% keepchrom),]
    dat.filt = dat.filt %>% ungroup()%>%
        mutate(start=floor(start/(10*binwin))*10*binwin) %>%
        group_by(chrom,start,lab)%>%
        summarize(count=sum(count))
    labs=unique(dat$lab)
    for (lab in labs){
        medcnt=median(dat.filt$count[which(dat.filt$lab==lab)])
        dat.filt$count[which(dat.filt$lab == lab &
                             dat.filt$count > 6 * medcnt )] = 6 * medcnt
        dat.filt$count[which(dat.filt$lab == lab)] =
            dat.filt$count[which(dat.filt$lab == lab)]/(medcnt/2)
    }
    dat.filt
}

output="/Users/kevin/Epigenomics/ambic-epigenome/data/day_cnv.pdf"
datpaths=c("/Users/kevin/Epigenomics/ambic-epigenome/data/HostD0.counts.10000.txt","/Users/kevin/Epigenomics/ambic-epigenome/data/StableGlutD0.counts.10000.txt","/Users/kevin/Epigenomics/ambic-epigenome/data/UnstableGlutD0.counts.10000.txt")

binwin=10000

faidxpath="/Users/kevin/Epigenomics/ambic-epigenome/data/picr_IgG2.fa.fai"
binwin=10000
processCNV <- function(datpaths,output,binwin,faidxpath){
    labels=sapply(strsplit(basename(datpaths),"[.]"),"[[",1)
    faidx=read_tsv(faidxpath,col_names=c("chrom","len","offset"))
    faidx$offset=ceiling(faidx$offset/binwin)*binwin
    dat.list=lapply(seq_along(datpaths),function(i){
        readcounts(datpaths[i])%>%
            mutate(lab=labels[i])})
    dat=do.call(rbind,dat.list)
    mtidx=grep("MT",dat$chrom)
    if (length(mtidx) > 0)dat=dat[-grep("MT",dat$chrom),] # remove MT chrom
    dat.bin=dat%>%
        mutate(start=floor(start/binwin)*binwin,
               end=ceiling(end/binwin)*binwin)
    dat.plt=dat.bin%>%select(-perbp)%>%
        spread(key=lab,value=count)%>%
        na.omit()
    # for genome-wide coverage plot
    dat.trim.plt=trimCov(dat.bin,faidx,binwin)
    dat.trim.plt$coordbin=dat.trim.plt$start+
        faidx$offset[match(dat.trim.plt$chrom,faidx$chrom)]
    dat.gather=dat.trim.plt %>% ungroup() %>%
        select(-chrom,-start)
    faidx.inc=faidx[which(faidx$chrom %in% dat.trim.plt$chrom),]
    seg.list=list()
    # segmentation?
    for ( lab in unique(dat.trim.plt$lab)){
        dat.sub=dat.trim.plt[which(dat.trim.plt$lab==lab),]
        cna=CNA(log(dat.sub$count),dat.sub$chrom,dat.sub$start,sampleid=lab)
        cna.sm=smooth.CNA(cna)
        cna.seg=segment(cna.sm,min.width=5,verbose=1)
        seg=as.tibble(segments.summary(cna.seg)[,c("chrom","loc.start","loc.end","seg.mean")])
        seg$seg.mean=exp(seg$seg.mean)
        seg$lab=lab
        seg.list[[lab]]=seg
    }
    dat.seg=do.call(rbind,seg.list)
    dat.seg$start=dat.seg$loc.start+faidx$offset[match(dat.seg$chrom,faidx$chrom)]
    dat.seg$end=dat.seg$loc.end+faidx$offset[match(dat.seg$chrom,faidx$chrom)]
    # rectangular plot
    multbin=1000*binwin
    dat.rect=dat.gather %>%
        mutate(bin=floor(coordbin/multbin)*multbin)%>%
        group_by(lab,bin)%>%
        summarize(count=mean(count))

    pdf(output,height=4,width=5,useDingbats=FALSE)
    combos=t(combn(labels,2))
    for (i in seq(dim(combos)[1])){
        onelab=combos[i,1]
        twolab=combos[i,2]
        plt.tb=dat.plt%>%
            select(one=onelab,two=twolab)
        g=ggplot(plt.tb,aes(x=one,y=two))+theme_bw()+
            geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
            labs(x=onelab,y=twolab)
        print(g)
    }
    for ( i in seq_along(datpaths)){
        samp=labels[i]
        plt.tb=dat.gather[which(dat.gather$lab==samp),]
        g = ggplot(plt.tb,aes(x=coordbin,y=count))+
            geom_point(size=0.2,alpha=0.2)+
            geom_vline(xintercept=faidx.inc$offset,
                       linetype="dashed",size=0.2,alpha=0.5)+
            theme_bw()+ ggtitle(samp)

        print(g)
    }
    # print(dim(dat.rect),count,lab)
    # print(bin,bin+multbin)
    g = ggplot(dat.rect,aes(ymin=0,ymax=count,
                      xmin=bin,xmax=bin+multbin,
                      fill=lab,group=lab))+
        geom_rect(alpha=0.4,colour=NA)+
        geom_vline(xintercept=faidx.inc$offset,
                   linetype="dashed",size=0.2,alpha=0.5)+
        theme_bw()+theme(legend.position="bottom")
    print(g)
    # g = ggplot(dat.seg)+
    #     geom_rect(aes(ymin=0,ymax=seg.mean,
    #                   xmin=start,xmax=end,fill=lab,group=lab),
    #               alpha=0.6,color=NA)+
    #     geom_vline(xintercept=faidx.inc$offset,
    #                linetype="dashed",size=0.2,alpha=0.5)+
    #     theme_bw()+theme(legend.position="bottom")
    # print(g)
    dev.off()
}



if (! interactive()){
    argsin=commandArgs(TRUE)
    args=parseArgs(argsin)
    print(args)
    processCNV(args$args,
               args$options$output,
               args$options$window,
               args$options$faidx)
}
