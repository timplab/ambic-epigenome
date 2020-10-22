#!/usr/bin/env python3
import sys
import os

in_fp = sys.argv[1]
width = 1e3 # 1kb for idt

# function to get overlap in intervals
# https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# since flanking is on 3' end, 
# if alignment strand is + : plasmid is inserted in 5' direction from the aligned region
# if alignment strand is - : plasmid is inserted in 3' direction from the aligned region
# so for sequencing into the plasmid from the genome, I need to design gRNA that :
# go toward 5' direction if alignment strand is + and
# go toward 3' direction if alignment strand is - and
# basically flip the alignment strand for gRNA strand
bp_list = []
with open(in_fp,"r") as fh :
    for line in fh :
        fields = line.strip().split("\t")
        contig = fields[3]
        strand = fields[4]
        if strand == "+" :
            start = int(fields[5])
            end = int(start + width - 1)
        elif strand == "-" :
            end = int(fields[6])
            start = int(end - width + 1)
        bp_list.append([contig,strand,start,end])

# iterate through each bp to group them
out_list = [bp_list[0]]
flag = 0
for bp in bp_list :
    for out in out_list :
        if bp[0] == out[0] and getOverlap(bp[2:],out[2:]) > 0 :
            flag = 1
    if flag == 0 : 
        out_list.append(bp)
    flag = 0

print("contig\talignment_strand\tstart\tend\tgRNA_direction\tcoord")
for out in out_list :
    if out[1] == "+" :
        direction = "left"
    elif out[1] == "-" :
        direction = "right"
    coord="{}:{}-{}".format(out[0],out[2],out[3])
    print("\t".join([str(x) for x in out] + [direction,coord]))
    

