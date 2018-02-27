#!/bin/bash

## generate bed files from bam files
## currently only SE option

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
    -g | --group) # group
      group=$2
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
bams=`readlink -f ${indir}/*${group}*nodup.bam`

##convert bam to tags
for bam in $bams;do
  s=${bam##*/}
  s=${s%.*}
  echo $s
  ## print as bed and shift
  bedtools bamtobed -i $bam | \
    awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";\
    if ($6 == "+"){$2=$2+4}\
    else if ($6 == "-"){$3=$3-5}\
    print $0}' | \
    sort -k1,1 -k2,2n |\
    gzip -nc > ${outdir}/${s}.bed.gz
done

## cross-cor analysis
## not applicable for SE data

## pool the beds
pool=${outdir}/${group}_pool.bed.gz
if [ -e $pool ];then
  rm $pool
fi
zcat ${outdir}/*${group}*.bed.gz |\
  sort -k1,1 -k2,2n |\
  gzip > ${outdir}/${group}_pool.bed.gz
