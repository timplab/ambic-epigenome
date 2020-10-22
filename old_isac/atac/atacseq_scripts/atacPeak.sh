#!/bin/bash

## generate bed files from bam files
## currently only SE option

##pipeline
atacpipe=~/Code/atac_dnase_pipelines

###argument parsing
while :
do
  case "$1" in
    -o | --outdir) # directory for final output
      mkdir -p $2
      outdir=`readlink -f $2`
      shift 2
      ;;
    -s | --sample) # name of the output
      samp=$2
      shift 2
      ;;
    -i | --input) #input dir
      indir=$2
      shift 2
      ;;
    -g | --gsize) ## genome size for macs
      g=$2
      shift 2
      ;;
    -pe | --paired-end) # flag for paired-end alignments
      pe=1
      shift 1
      ;;
    *) break
      ;;
  esac
done
## input
bed=`readlink -f ${indir}/${samp}*bed.gz`
## default for ENCODE
pval=0.01
smoothwin=73
shiftn=$((($smoothwin+1)/2))
## output
outpre=${outdir}/${samp}

## limit the memory usage
#ulimit -v 30000000

echo "sample name : $samp"
macs2 callpeak \
  -t $bed -f BED -n $outpre -g $g -p $pval \
  --nomodel --shift $shiftn --extsize $smoothwin \
  -B --SPMR --keep-dup all --call-summits 2> ${outpre}_peakcalling.log

## Sort by Col8 in descending order and 
# replace long peak names in Column 4 with Peak_<peakRank>
sort -k 8gr,8gr ${outpre}_peaks.narrowPeak | \
  awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | \
  gzip -nc > ${outpre}.narrowPeak.gz
