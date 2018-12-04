#! /usr/bin/env python
import os
import sys
import argparse
import math

def parseArgs() :
    # parse args
    parser = argparse.ArgumentParser( description='bin genome into equal sized chunks' )
    parser.add_argument('-f', '--fai',type=argparse.FileType('r'),required=True,
            help="reference genome index")
    parser.add_argument('-b', '--binsize', type=int, required=False,
            default=10000, help="binning window size")
    parser.add_argument('--filetype',type=str,required=False,default="bed",
            help="output file type, bed or coord")
    args = parser.parse_args()
    return args

def output_bed(contig,window) :
    print("{}\t{}\t{}".format(contig,window[0],window[1]),file=sys.stdout)

def output_coord(contig,window) :
    print("{}:{}-{}".format(contig,window[0]+1,window[1]),file=sys.stdout)

def main() :
    args = parseArgs()
    if args.filetype == "bed" :
        outfunc = output_bed
    elif args.filetype == "coord" :
        outfunc = output_coord
    for line in args.fai :
        fields = line.strip().split()
        contig = fields[0]
        contigsize = int(fields[1])
        numbins = math.ceil(contigsize/args.binsize)
        wins = [ [x*args.binsize,(x+1)*args.binsize] for x in range(numbins) ]
        wins[-1][1] = contigsize 
        [ outfunc(contig,window) for window in wins ]
if __name__=="__main__":
    main()

