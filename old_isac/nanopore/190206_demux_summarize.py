#!/home/ubuntu/Code/miniconda3/bin/python
import os
import numpy as np

max_num = 100000
os.chdir("/shared/demux_summary")
sums = [ x for x in os.listdir() if "summary.txt" in x ]

n_list = list()
s_list = list()
for fp in sorted(sums) :
    sample = fp.split(".")[0]
    in_fh = open(fp,'r')
    name_list = list()
    num_list = list()
    for line in in_fh:
        fields = line.strip().split()
        name_list.append(fields[1])
        num_list.append(int(fields[0]))
    num_list = np.array(num_list)
    total = np.sum(num_list)
    order_idx = np.argsort(num_list)[::-1]
    name_ordered = np.array(name_list)[order_idx]
    num_ordered = num_list[order_idx]
    name_top = name_ordered[0:4]
    num_top = num_ordered[0:4]
    order_idx = np.argsort(name_top)
    name_top = list(name_top[order_idx])
    num_top = list(num_top[order_idx])
    n_list.append('\t'.join([str(x) for x in [sample,total]+num_top[0:3]]))
    s_list.append('\t'.join([str(x) for x in [sample,total]+name_top[0:3]]))

for x in n_list :
    print(x)
for x in s_list :
    print(x)

