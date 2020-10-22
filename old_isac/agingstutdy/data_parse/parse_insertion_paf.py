#!/usr/bin/env python3
import sys
import os
import numpy as np

in_fp = sys.argv[1]
align_fp = sys.argv[2]

# function to get overlap in intervals
# https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
def annotate_read(paf_list) :
    # making a list of lists
    # each row will have : qname, qlen, qstart, qend, rname, rstart, rend, direction
    anno_list = []
    qlen = int(paf_list[0][1])
    qname = paf_list[0][0]
    for paf in paf_list :
        qstart = int(paf[2])
        qend = int(paf[3])
        direction = paf[4]
        rname = paf[5]
        rstart = int(paf[7])
        rend = int(paf[8])
        anno_list.append([qname,qlen,qstart,qend,rname,rstart,rend,direction])
    return anno_list

# first parsing plasmid alignments
rlen_dict = dict()
reads_dict = dict()
hits_dict = dict()
paf_dict = dict()
# get coordinates and query matches
# let's also record subject hits
with open(in_fp,"r") as fh :
    for line in fh :
        fields = line.strip().split("\t")
        qname = fields[0]
        # add the raw data to paf_dict
        if qname not in paf_dict.keys() :
            paf_dict[qname] = []
        paf_dict[qname].append(fields)
        # get query coords and add to reads_dict
        query_coord = [ int(fields[2]), int(fields[3]) ]
        if qname not in reads_dict.keys() :
            reads_dict[qname] = []
        reads_dict[qname].append(query_coord)
        rlen_dict[qname] = int(fields[1])
        # subject hits
        subject_coord = [ int(fields[7]),int(fields[8]) ]
        if qname not in hits_dict.keys() :
            hits_dict[qname] = []
        hits_dict[qname].append(subject_coord)
# let's get total width of alignments
matchlen_dict = dict()
for key in reads_dict.keys() :
    query_coords = reads_dict[key]
    pos_list = [list(range(x[0],x[1])) for x in query_coords]
    pos = set([ x for y in pos_list for x in y ])
    matchlen_dict[key] = len(pos)

# maximum of the total match in host cells is ~1.4kb
# so let's require at least 2kb total match to consider a real hit
qnames_hit = [ key for key in matchlen_dict.keys() if matchlen_dict[key] > 2e3 ]
if len(qnames_hit) == 0 : 
    print("No reads with valid hits", file = sys.stderr)
    sys.exit()
else :
    print("{} reads with valid hits".format(len(qnames_hit)),file = sys.stderr)

###########################################
# Start parsing valid hits
###########################################
# get min start and max end
minmax_dict = dict()
for key in qnames_hit :
    coords = reads_dict[key]
    minstart = min([x[0] for x in coords ])
    maxend = max([x[1] for x in coords ])
    minmax_dict[key] = [minstart,maxend]
    al_list = annotate_read(paf_dict[key])
    for al in al_list :
        print("\t".join([str(x) for x in al]))

###########################################
# Flanking reads
###########################################
# parse flanking - how much of non-plasmid sequence on each side of the query?
# we also don't want reads that are solely on plasmid
thr = 1e3
# each list will have [ left flanking bps, right flanking bps, length match ]
flank_dict = dict()
for key in minmax_dict.keys() :
    left = minmax_dict[key][0]
    end = minmax_dict[key][1]
    rlen = rlen_dict[key]
    right = rlen - end
    matchlen = matchlen_dict[key]
    if left < thr and right < thr :
        continue
    # if both sides are flanking, then unless at least one whole plasmid is in the match, this is error
    if left > thr and right > thr and matchlen < 10e3 :
        continue
    flank_dict[key] = [left,right, matchlen]
qnames = flank_dict.keys()

if len(qnames) == 0 :
    print("No flanking reads", file = sys.stderr)
    sys.exit()
else :
    print("{} flanking reads".format(len(qnames)), file = sys.stderr)

###########################################
# Parse flanking reads
###########################################
flankcoord_dict = { x:[] for x in flank_dict.keys() }
for key in flank_dict.keys() :
    rlen = rlen_dict[key]
    left = flank_dict[key][0]
    right = flank_dict[key][1]
    if left > thr :
        coord = [0,left]
        flankcoord_dict[key].append(coord)
    if right > thr :
        coord = [rlen - right,rlen]
        flankcoord_dict[key].append(coord)
      
# let's read in the reference alignment paf
ref_dict = dict()
with open(align_fp,"r") as fh :
    for line in fh :
        fields = line.strip().split("\t")
        qname = fields[0]
        if qname in flank_dict.keys() :
            qstart = int(fields[2])
            qend = int(fields[3])
            contig = fields[5]
            rstart = int(fields[7])
            rend = int(fields[8])
            if qname not in ref_dict.keys() :
                ref_dict[qname] = []
            ref_dict[qname].append([qstart,qend,fields[5],fields[4],rstart,rend])

# get reference coordinates
for key in ref_dict.keys() :
    flanks = flankcoord_dict[key]
    for align in ref_dict[key]:
        qcoords = [align[0],align[1]]
        for flank in flanks :
            if getOverlap(qcoords,flank) > 0 :
                print("\t".join([key] +[str(x) for x in align]),file = sys.stderr)

