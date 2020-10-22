#!/usr/bin/env python3
import sys

fqnumfp = sys.argv[1]
sumnumfp = sys.argv[2]

# one
def parse_numfile(fp) :
    num_dict = dict()
    with open(fp,"r") as fh :
        for line in fh :
            fields = line.strip().split(" ")
            num = int(fields[0])
            sample = fields[1].split("/")[-1].split(".")[0]
            num_dict[sample] = num
    return num_dict

fqnum_dict = parse_numfile(fqnumfp)
sumnum_dict = parse_numfile(sumnumfp)

for key in fqnum_dict.keys() :
    sumnum = sumnum_dict[key]-1
    fqnum = fqnum_dict[key]
    if sumnum != fqnum :
        print("Error in {}".format(key))

