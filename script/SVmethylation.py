#! /usr/bin/env python
import os
import math
import sys
import argparse
import gzip
import numpy as np
from methylbed_utils import MethRead,SnifflesEntry,make_coord,read_bam
import pysam
import re
import multiprocessing as mp
from multiprocess_utils import listener,init_mp,close_mp
import time
start_time = time.time()

def parseArgs() :
    # dir of source code
    srcpath=sys.argv[0]
    srcdir=os.path.dirname(os.path.abspath(srcpath))
    # parser
    parser = argparse.ArgumentParser(description='parse methylation around SVs')
    parser.add_argument('-t','--threads',type=int,required=False,default=2, 
            help="number of parallel processes (default : 2 )")
    parser.add_argument('-v','--verbose', action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-s','--sniffles',type=argparse.FileType('r'),required=False, 
            default=sys.stdin,help="sniffles vcf output (default stdin)")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="bam file")
    parser.add_argument('-c','--cpg',type=os.path.abspath,required=True,
            help="gpc methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('-w','--window',type=int,required=False,
            default=200,help="window for methylation")
    parser.add_argument('-o','--output',type=argparse.FileType('w'),required=False, 
            default = sys.stdout,help="output path (default : stdout)")
    # parse args
    args = parser.parse_args()
    args.srcdir=srcdir
    return args

def read_tabix(fpath,window) :
    with pysam.TabixFile(fpath) as tabix :
        entries = [x for x in tabix.fetch(window)]
    reads = [MethRead(x) for x in entries]
    rdict = dict()
    for read in reads :
        try :
            rdict[read.qname].append(read)
        except :
            rdict[read.qname] = [read]
    return rdict

def getRegMeth(read_list,start,end) :
    callarray = np.concatenate([ x.callarray for x in read_list])
    regidx = np.where(np.logical_and(callarray[:,0]>=start, callarray[:,0]<=end))
    callreg = callarray[regidx,1].flatten()
    sigidx = np.argwhere(callreg != -1)
    sigreg = callreg[sigidx]
    methcount = np.count_nonzero(sigreg)
    return len(sigreg),methcount
    
def parse_methylation(q,sv,cpg,gpc,start,end,tag) :
    qname = cpg[0].qname
    taglist = tag.split("_")
    cpgcov,cpgmeth = getRegMeth(cpg,start,end)
    gpccov,gpcmeth = getRegMeth(gpc,start,end)
    if (gpccov == 0 or cpgcov == 0) : return
    line = '\t'.join([str(x) for x in [ sv.chrom,sv.pos,sv.pos,
        sv.info["CHR2"],sv.info["END"],sv.info["END"],
        qname,sv.id,".",".",taglist[1],taglist[0],cpgmeth,gpcmeth,cpgcov,gpccov]])
    q.put((qname+sv.id+tag,line))

def TRA_methylation(sv,bamfn,cpgfn,gpcfn,methwin,verbose,q) :
    # windows for fetching reads - 
    win = 300
    win1 = make_coord(sv.chrom,sv.pos-win,sv.pos+win)
    win2 = make_coord(sv.info["CHR2"],sv.info["END"]-win,sv.info["END"]+win)
    # fetch reads
    try : 
        bam_dicts = [ read_bam(bamfn,w) for w in [win1,win2] ]
    except ValueError : 
        return
    cpg_dicts = [ read_tabix(cpgfn,w) for w in [win1,win2] ]
    gpc_dicts = [ read_tabix(gpcfn,w) for w in [win1,win2] ]
    qnames = list(bam_dicts[0].keys()) + list(bam_dicts[1].keys())
    for qname in set(qnames) :
        if (qname in sv.rnames or 
                ( qname in bam_dicts[0] 
                    and qname in bam_dicts[1] )):
            # this read is an SV and has both parts
            if qname in cpg_dicts[0].keys() and qname in gpc_dicts[0].keys() :
                cpg1 = cpg_dicts[0][qname]
                gpc1 = gpc_dicts[0][qname]
                parse_methylation(q,sv,cpg1,gpc1,sv.pos-methwin,sv.pos+methwin,"destination_SV")
            if qname in cpg_dicts[1].keys() and qname in gpc_dicts[1].keys() :
                cpg = cpg_dicts[1][qname]
                gpc = gpc_dicts[1][qname]
                start,end,tag = (sv.info["END"]-methwin,sv.info["END"]+methwin,"origin_SV")
            else :
                continue
        elif ( qname in bam_dicts[1] ) :
            # non-SV origin
            if qname in cpg_dicts[1].keys() and qname in gpc_dicts[1].keys() :
                cpg = cpg_dicts[1][qname]
                gpc = gpc_dicts[1][qname]
            else : 
                continue
            coords = [ pos for x in bam_dicts[1][qname] for pos in [x.reference_start,x.reference_end] ]
            start,end = (sv.info["END"]-methwin,sv.info["END"]+methwin)
            tag = "origin_nonSV"
        elif ( qname in bam_dicts[0] ) :
            # non-SV destination
            try :
                cpg = cpg_dicts[0][qname]
                gpc = gpc_dicts[0][qname]
            except :
                continue
            coords = [ pos for x in bam_dicts[0][qname] for pos in [x.reference_start,x.reference_end] ]
            start,end = (sv.pos-methwin,sv.pos+methwin)
            bp = [ x for x in coords if x >= sv.pos-win and x <= sv.pos+win ]
            tag = "destination_nonSV"
        parse_methylation(q,sv,cpg,gpc,start,end,tag)

if __name__=="__main__":
    args=parseArgs()
    svlines = [SnifflesEntry(x) for x in args.sniffles.readlines() if x[0]!="#"]
    svreads = [ read for x in svlines for read in x.rnames]
    fetchwin = 2000
    if len(svreads) == 0 : 
        windows = [ "scaffold_16:15300006-15310006" ]
    else : 
        windows = [ make_coord(x.chrom,x.pos-fetchwin,x.pos+fetchwin) for x in svlines ]
    meth = read_tabix(args.cpg,windows[0])
#    bams = read_bam(args.bam,windows[0])
    svmeth = list()
    if len(svreads) >= 0 :  # change this to >= to make it inclusive
        for qname in meth.keys() :
            svmeth.append(meth[qname])
    else : 
        for qname in svreads :
            if qname in meth.keys() :
                svmeth.append(meth[qname])
    for readlist in svmeth :
        for read in readlist :
            print("\t".join(read.fields),file=args.output)
                    
    if args.verbose : print("time elapsed : {} seconds".format(time.time()-start_time),file=sys.stderr)
