#!/usr/bin/Rscript
# plot average curves
# input is a tsv file that has sample,replicate,readnum,total-bases,N50,Q1,mean,Q3,max
library(optparse)
library(tidyverse)

parseArgs = function(argsin){
    arglist=list(
        make_option(c("-o","--output"),type="character",default=NULL,
                    help="output prefix",metavar="/path/to/out")
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

readsum <- function(fpath){
    cnames=c("filename","readname","lane","channel","qual","idk","len","idk2","rlen")
    dat=read_tsv(fpath,col_names=cnames)
    dat
}

out="/dilithium/Data/Nanopore/Analysis/180823_qc_lambda/180823_qc_lambda_QC.pdf"
datpath="/dilithium/Data/Nanopore/oxford/180823_qc_lambda/180823_qc_lambda.tsv.gz"
ontDistro <- function(datpath,out){
    dat=readsum(datpath)
    dat.filt = dat %>%
        filter(rlen<60000)
    # plot
    g = ggplot(dat,aes(x=rlen))+
        theme_bw()+
        scale_y_log10()+
        expand_limits(x=0)
    g.hist = g + geom_histogram()
    g.freqpoly = g + geom_freqpoly()
    g.freq.log = g.freqpoly + scale_x_log10()
    g.hist.log = g.hist + scale_x_log10()
    gd = ggplot(dat, aes(x=rlen,y=..scaled..))+
        scale_y_log10()+
        geom_density(alpha=0.5) + theme_bw()
    gd.log = gd + scale_x_log10()
        
    # output into a pdf file
    pdf(out,height=4,width=6,useDingbats=FALSE)
    print(g.freq.log)
    print(g.hist.log)
    dev.off()

}

if (! interactive()){
    argsin=commandArgs(TRUE)
    args=parseArgs(argsin)
    print(args)
    ontDistro(args$args,args$options$output)
}
