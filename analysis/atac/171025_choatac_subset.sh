#!/bin/bash

script=/home/isac/Code/ilee/atac/atacseqPipeline.py
wdir=/dilithium/Data/NGS/Aligned/171025_choatac/subset
btidx=/mithril/Data/NGS/Reference/cho/criGri1_plasmid2
kdir=${wdir}/choatac_subset_metatdata.csv

$script -k $kdir -o $wdir -x $btidx -g mm -m peak
