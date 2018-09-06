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

readqc <- function(fpath){
    cnames=c("samp","rep","read_num","total_bp","n50","q1","mean_length","q3","max_length")
    dat=read_tsv(fpath,col_names=cnames)
    dat
}

outpre="/home/isac/Dropbox/Data/ambic/qc/QClambda.pdf"
datpath="/home/isac/Dropbox/Data/ambic/qc/ontqc.txt"
ontQC <- function(datpath,outpre){
    dat=readqc(datpath)
    print(dat)
    # get numbers for simulating yields
    yrange=c(min(dat$mean_length),max(dat$mean_length)) # take min and max of mean read lengths
    yields=round(seq(ceiling(min(dat$total_bp)/1e9),
               floor(max(dat$total_bp)/1e9),
               length.out=5)) # get 5 yield numbers based on min and max of yields from runs
    n=100 # make 100 data points per simulated yield
    rlens=seq(yrange[1],yrange[2],length.out=n)
    sim.tb=tibble(rlen=rlens)
    for ( yield in yields ){
        lab=as.character(yield)
        sim.tb[,lab]=yield*1e9/sim.tb$rlen
    }
    # max matrix of values into vertical array
    sim.spread=sim.tb %>%
        gather(yield,rnum,-rlen) %>%
        mutate(yiled=as.numeric(yield)) %>%
        filter(rnum<=max(dat$read_num)) # remove huge read num points that exceed max from data
    # plot
    g1 = ggplot(dat, aes(read_num, mean_length))+
        geom_point(aes(color=samp,shape=factor(rep)),size = 2) +
        theme_bw() + xlab("Number of Reads") + 
        ylab("Average Read Length") +
        ggtitle("Results of Sequencing Runs")
    g2 = g1 + geom_line(data = sim.spread,
                        mapping = aes(rnum,rlen,group = yield),
                        color = "orange",size = 0.5,
                        linetype = "dashed") +
        guides(color=guide_legend(paste(yields,collapse="\t")))

    # output into a pdf file
    pdf(paste0(outpre,"_yields.pdf"),height=4,width=6,useDingbats=FALSE)
    print(g2)
    dev.off()

}

if (! interactive()){
    argsin=commandArgs(TRUE)
    args=parseArgs(argsin)
    print(args)
    ontQC(args$args,args$options$output)
}
