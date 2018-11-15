#!/bin/bash
root=/shared
srcpath=$root/Code/ilee/oxford
slurmpath=$srcpath/slurm
rawdir=$root/raw
callroot=${root}/called
[ -e $callroot ]||mkdir $callroot
fqdir=${root}/fqchunks

for tag in host IgG
do
  for rep in 3 2 1
  do
    fast5path=`ls -d $rawdir/*$tag*$rep*`
    base=$(basename "$fast5path")
    calldir=$callroot/$base
    [ -e $calldir ]||mkdir -p $calldir
#    ## identify fast5 folder path and perform basecalling
#    echo "basecalling reads in $fast5path"
#  echo "$slurmpath/call_wrapper.sh -i $fast5path -o $calldir \
#      -f FLO-MIN106 -s $slurmpath &> $callroot/callwrapper_$base.log"
#  $slurmpath/call_wrapper.sh -i $fast5path -o $calldir \
#      -f FLO-MIN106 -s $slurmpath &> $callroot/callwrapper_$base.log
  # cat the files  
  for dir in `find $calldir/* -maxdepth 0 -type d`;do
    ind=$(basename "$dir")
    prefix=${fqdir}/$base.$ind
    echo "concatentating $base index $ind fastq files"
    logpath=${prefix}.cat.log
    sbatch -e $logpath -o $logpath $slurmpath/cat_bcall.scr -i $dir -p $prefix
  done
  done
done



