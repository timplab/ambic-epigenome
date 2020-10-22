#!/bin/bash
export LC_ALL=C
root=/mnt/c/Data
TMP=$root/tmp
TEMP=$root/tmp
TMPDIR=$root/tmp
export TMP TEMP TMPDIR
mdir=$root/mcall
bdir=$root/mbed
pooldir=$root/pooled_rep/mbed
day=0
scr=/home/ubuntu/Code/ilee/nanopolish/mtsv2bedGraph.py

if [ "$1" == "mbed" ];then
  for tsv in $(find $mdir -name "*tsv.gz"); do
    base=$(basename "$tsv")
    base=${base%%.*}
    bases="$bases $base"
    mkdir $root/tmp/$base
  done
  mtsv="$mdir/{}.cpg.meth.tsv.gz"
  mbed="$bdir/{}.cpg.meth.bed.gz"
  com="gunzip -c $mtsv | python $scr |\
    sort -T $root/tmp/{} -k1,1 -k2,2n | bgzip > $mbed"
  parallel $com ::: $bases
fi

if [ "$1" == "pool" ];then
  for cell in Host StableGln StableNogln UnstableGln UnstableNogln; do
    for rep in 1 2 3; do
      samp="CHOZN${cell}Day${day}_$rep"
      samps="$samps $samp"
      mkdir $root/tmp/$samp
    done
  done
  samp="{}"
  mbeds="$bdir/{}*cpg.meth.bed.gz"
  out="$pooldir/{}.cpg.meth.bed.gz"
  com="gunzip -c $mbeds | sort -T $root/tmp/{} -k1,1 -k2,2n |\
    bgzip > $out"
#  com="echo $mbeds $out"
  parallel $com ::: $samps
fi

if [ "$1" == "getfreq" ];then
  scr=/home/ubuntu/Code/nanoNOMe/scripts/parseMethylbed.py
  for cell in Host StableGln StableNogln UnstableGln UnstableNogln; do
    for rep in 1 2 3; do
      samp="CHOZN${cell}Day${day}_$rep"
      samps="$samps $samp"
      mkdir $root/tmp/$samp
    done
  done
  samp="{}"
  mbed="$pooldir/{}.cpg.meth.bed.gz"
  out="$pooldir/../mfreq/{}.cpg.meth.freq.txt.gz"
  com="python $scr frequency -i $mbed |\
    bgzip > $out && tabix -b 2 -e 2 $out"
  parallel $com ::: $samps
fi
