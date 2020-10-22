#!/usr/bin/env python3
# halt
import os
import sys
import argparse
import pysam
from collections import namedtuple
import re
import numpy as np

def parseArgs() :
    parser = argparse.ArgumentParser(description='locate the site of insertion')
    parser.add_argument('-v','--verbose',action='store_true',default=False,
            help="verbose output")
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help='sorted and indexed bam file')
    parser.add_argument('-c','--contig',nargs='*',required=True,
            help="name of the insertion contig")
    parser.add_argument('--bin',type=int,required=False,default=200,
            help="binning width (default 200)")
    args = parser.parse_args()
    return args

def bool_to_strand(logical) :
    if logical :
        return "-"
    else :
        return "+"

def get_query_start(cigar) :
    try : 
        return int(re.split('H',cigar)[0])
    except ValueError :
        return 0 

alignMetric = namedtuple('alignMetric',['rname','rstart','rend','strand','qstart','qend','supp'])
alignSummary = namedtuple('alignSummary',['qstart','qend','rname','rstart','rend','strand'])
if __name__=="__main__":
    args = parseArgs()
    def numbin(number):
        return int(number/args.bin)*args.bin
    bam = pysam.AlignmentFile(args.bam,'rb')
    # first get names of the query
    qlist = list()
    rlen_list = list()
    for contig in args.contig :
        for read in bam.fetch(contig) :
            if not read.has_tag("SA") :
                continue
            qlist.append(read.query_name)
            rlen_list.append(read.infer_read_length())
    queryname = qlist[np.argmax(rlen_list)]
    # cycle through alignments
    qset = set(qlist)
    adict = dict()
    for read in bam.fetch() :
        qname = read.query_name 
        if qname in qset :
            if qname not in adict.keys() :
                adict[qname] = [read]
            else :
                adict[qname].append(read)
    # cycle through reads
#    alist = list(adict.values())
#    amax = alist[np.argmax([ len(x) for x in alist ])]
#    aq = amax[0].query_name
    for qname in adict.keys():
        summary_list = list()
        read0 = None
        alignments = adict[qname]
        if len(alignments) == 1 : continue # skip reads with only one alignment
        ametrics = list()
        # get the qstart and order along query
        for read in alignments :
            rlen = read.infer_read_length()
            offset = get_query_start(read.cigarstring) 
            qstart = offset + read.query_alignment_start
            qend = offset + read.query_alignment_end
            ametrics.append(alignMetric(read.reference_name,read.reference_start,
                    read.reference_end,bool_to_strand(read.is_reverse),
                    qstart,qend,read.is_supplementary))
        # parse split reads
        qstart_list = [ x.qstart for x in ametrics ]
        order_idx = sorted(range(len(qstart_list)),key=qstart_list.__getitem__)
        for i in order_idx :
            if read0 is None :
                read0 = ametrics[i]
                continue
            read1 = ametrics[i]
            range0 = set(range(read0.qstart,read0.qend))
            range1 = set(range(read1.qstart,read1.qend))
            ovl = len(range0 & range1)
            avgrange = (len(range0)+len(range1))/2
            # same part of the query aligned 
            if ovl/avgrange > 0.5 :
                if len(range0) > len(range1) :
                    rl = [ read0, read1 ]
                else : 
                    rl = [ read1, read0 ]
                try : read1 = [ read for read in rl if read.supp is False ][0]
                except IndexError :
                    read1 = rl[0]
            else :
                summary_list.append(alignSummary(read1.qstart,read1.qend,
                        read1.rname,read1.rstart,read1.rend,read1.strand))
            read0 = read1
        if len(summary_list) <= 1 : continue
        for x in summary_list :
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                qname,x.qstart,x.qend,x.rname,x.rstart,x.rend,x.strand,rlen))

    bam.close()
