#!/usr/bin/env python3
import sys

qcpath = sys.argv[1]

# all the samples
cells = [ "Host","StableGln","StableNogln","UnstableGln","UnstableNogln" ]
days = [ 0,30,60,90 ]
reps = [ 1,2,3 ]
samples = list()
for cell in cells :
    for day in days :
        for rep in reps :
            sample = "CHOZN{}Day{}_{}".format(cell,day,rep)
            samples.append(sample)

qc_dict = {key:[0,0] for key in samples}
with open(qcpath,'r') as fh :
    fh.readline()
    for line in fh :
        fields = line.strip().split(',')
        sample = '_'.join(fields[0].split('/')[-1].split('_')[0:2])
        qc_dict[sample][0] += int(fields[1])
        qc_dict[sample][1] += int(fields[2])
    print("sample,numreads,yieldGb,avglen")
    for key in sorted(qc_dict.keys()) :
        nums = qc_dict[key]
        numbase_gb = int(nums[1]/1e9)
        avglen = int(round(nums[1]/(nums[0]+1e-9)))
        print("{},{},{},{}".format(
            key,nums[0],numbase_gb,avglen))
