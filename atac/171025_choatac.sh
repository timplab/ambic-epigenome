#!/bin/bash

script=/home/isac/Code/ilee/atac/atacseqPipeline.py
wdir=/dilithium/Data/NGS/Aligned/171025_choatac/
btidx=/mithril/Data/NGS/Reference/cho/picr_IgG2/picr_IgG2
#kdir=${wdir}/choatac_metatdata.csv
kdir=${wdir}/choatac.tmp.csv

$script -k $kdir -o $wdir -x $btidx -g mm -m bamprocess
