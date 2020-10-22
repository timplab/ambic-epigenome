#!/usr/bin/env python3
import sys
import os

in_fp = sys.argv[1] # plasmid alignments paf
align_fp = sys.argv[2] # genome alignments paf

rlen_dict = dict()
reads_dict = dict()
# get coordinates and query matches
with open(in_fp,"r") as fh :
    for line in fh :
        fields = line.strip().split("\t")
        qname = fields[0]
        query_coord = [ int(fields[2]), int(fields[3]) ]
        if qname not in reads_dict.keys() :
            reads_dict[qname] = []
        reads_dict[qname].append(query_coord)
        rlen_dict[qname] = int(fields[1])

# get min start and max end
minmax_dict = dict()
for key in reads_dict.keys() :
    coords = reads_dict[key]
    minstart = min([x[0] for x in coords ])
    maxend = max([x[1] for x in coords ])
    minmax_dict[key] = [minstart,maxend]

# which ones have enough flanking?
# since adaptors bind to the cas9 cut site, all of our flanking regions should be on the 3' end
# but let's verify it / be generous with considering flanking sites
# in running this with our first round data, I have confirmed that indeed, all flanking sites are on the 3' end
thr = 300
flank_dict = dict()
for key in minmax_dict.keys() :
    left = minmax_dict[key][0]
    end = minmax_dict[key][1]
    rlen = rlen_dict[key]
    right = rlen - end
    # In whole genome sequencing, both sides can have genomic seq, so still use this
    if left > thr :
        flank_dict[key] = [[0,left]]
    if right > thr :
        if key in flank_dict.keys() :
            flank_dict[key].append([end-1,rlen])
        else : 
            flank_dict[key] = [[end-1,rlen]]

# let's read in the reference alignment paf
ref_dict = dict()
with open(align_fp,"r") as fh :
    for line in fh :
        fields = line.strip().split("\t")
        qname = fields[0]
        if qname in flank_dict.keys() :
            qstart = int(fields[2])
            qend = int(fields[3])
            if qend - qstart < thr : continue
            contig = fields[5]
            rstart = int(fields[7])
            rend = int(fields[8])
            if qname not in ref_dict.keys() :
                ref_dict[qname] = []
            ref_dict[qname].append([qstart,qend,fields[5],fields[4],rstart,rend])

# function to get overlap in intervals
# https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# compare
# since flanking is on 3' end, 
# if alignment strand is + : plasmid is inserted in 5' direction from the aligned region
# if alignment strand is - : plasmid is inserted in 3' direction from the aligned region
# so for sequencing into the plasmid from the genome, I need to design gRNA that :
# go toward 5' direction if alignment strand is + and
# go toward 3' direction if alignment strand is - and
# basically flip the alignment strand for gRNA strand
for key in ref_dict.keys() :
    flanks = flank_dict[key]
    for align in ref_dict[key]:
        qcoords = [align[0],align[1]]
        for flank in flanks :
            if getOverlap(qcoords,flank) > 0 :
                print("\t".join([key] +[str(x) for x in align]))

