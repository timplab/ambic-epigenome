#!/bin/bash

## generate bed files from bam files
## currently only SE option

###argument parsing
while :
do
  case "$1" in
    -i | --in)
      in=$2
      shift 2
      ;;
    -o | --out) # directory for final output
      out=$2
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
##convert bam to tags
bedtools bamtobed -i $in| \
  awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
  sort -k1,1 -k2,2n |\
  bgzip > $out

