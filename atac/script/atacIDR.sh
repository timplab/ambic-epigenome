#!/bin/bash

## IDR
## unfinished

##pipeline
atacpipe=~/Code/atac_dnase_pipelines

###argument parsing
while :
do
  case "$1" in
    -i | --indir)
      indir=$2
      shift 2
      ;;
    -o | --outdir) # directory for final output
      outdir=`readlink -f $2`
      shift 2
      ;;
    -b | --base) # bases in quotation
      base=$2
      shift 2
      ;; 
    *) break
      ;;
  esac
done
## input

echo "$0"
for samp in $base
do
  peak=`readlink -f $indir/${samp}.narrowPeak.gz`

done
