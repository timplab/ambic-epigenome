#! /usr/bin/env python

## Create a binned coverage matrix using a bed file. 

import math
import sys
import csv
import argparse
import gzip
from collections import namedtuple

def parseArgs():
    parser = argparse.ArgumentParser( description='bin read counts from bed file')
    parser.add_argument('-i', '--input', type=str, required=False,
            help="input coverage bedGraph file")
    parser.add_argument('-w','--binwidth',type=int,required=False,default=100000,
            help="window for binning")
    args = parser.parse_args()
    return args

class bedQuery:
    def __init__(self,string):
        self.fields=string.decode().strip().split("\t")
        self.rname=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.qname=self.fields[3]
        self.mapq=int(self.fields[4])
        self.strand=self.fields[5]
class covWindow:
    def __init__(self,query,width):
        self.rname=query.rname
        self.start=query.start
        self.width=width
        self.end=query.end
        self.binend=self.start+self.width
        self.cov=1
        self.perbp=0
    def update(self,newquery):
        if ( newquery.start < self.binend and
                newquery.rname == self.rname ):
            self.cov+=1
            return
        else :
            self.perbp = self.cov / (self.end-self.start)
            self.printCov()
            self.__init__(newquery,self.width)
        return
    def printCov(self):
        print("\t".join([str(x) for x in [self.rname,
            self.start,
            self.end,
            self.cov,
            self.perbp]]))


def binCoverage(in_fp,width=25000):
    if in_fp:
        in_fh = gzip.open(in_fp,'r')
    else:
        in_fh = sys.stdin
    for line in in_fh:
        query=bedQuery(line)
        try :
            window.update(query)
        except NameError :
            window=covWindow(query,width)
    in_fh.close()

if __name__=="__main__":
    args=parseArgs()
    binCoverage(args.input,args.binwidth)
