#!/usr/bin/python
import sys
import os
summary_path = sys.argv[1]
print("barcode\tread_number\ttotal_length")
reads_dict = dict()
with open(summary_path,'r') as fh :
    fh.readline()
    for line in fh :
        fields = line.strip().split("\t")
        barcode = fields[19]
        try : 
            reads_dict[barcode].append(int(fields[12]))
        except KeyError : 
            reads_dict[barcode] = [int(fields[12])]
    for key in reads_dict.keys() :
        l = reads_dict[key]
        n = len(l)
        s = sum(l)
        print("\t".join([key,str(n),str(s)]))
