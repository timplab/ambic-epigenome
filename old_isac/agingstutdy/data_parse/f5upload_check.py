#!/home/ubuntu/Code/miniconda3/bin/python
import sys
import os
import math

fq_fh = open("/mnt/c/qc/190320_fastq_qc.csv","r")
f5_fh = open("/mnt/c/qc/numf5.txt","r")

fq_dict = dict()
fq_fh.readline()
for line in fq_fh : 
    fields = line.strip().split(",")
    fqnum = int(fields[1])
    f5num = math.ceil(fqnum/4000)
    fq_dict[fields[0]] = f5num

f5_dict = dict()
for line in f5_fh :
    fields = line.strip().split(",")
    f5_dict[fields[0]] = int(fields[1])

for key in f5_dict.keys() :
    f5num = f5_dict[key]
    predictnum = fq_dict[key]
    print("{} predicted, {} present, {} more than predicted".format(
        predictnum,f5num,f5num-predictnum))
    if f5num-predictnum < 0 : print("warning")


f5_fh.close()
fq_fh.close()
