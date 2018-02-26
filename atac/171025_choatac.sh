#!/bin/bash

script=./script/atacseqPipeline.py
wdir=/dilithium/Data/NGS/Aligned/171025_choatac/
btidx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd
kdir=${wdir}/choatac.metadata.csv
kdir=${wdir}/one.csv

if [ "$1" == "trim" ];then
  $script -k $kdir -o $wdir -x $btidx -g mm -m trim
fi
if [ "$1" == "align" ];then
  $script -k $kdir -o $wdir -x $btidx -g mm -m align
fi
