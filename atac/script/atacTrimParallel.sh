#!/bin/bash

## Usage: ./atacTrim.sh -o [outdir] \
##  -pe (if PE) -i [inputdir] -s [sampleID]

###argument parsing
while :
do
  case "$1" in
    -o | --outdir) # directory for final output
      mkdir -p $2
      outdir=`readlink -f $2`
      shift 2
      ;;
    -pe | --paired-end) #is it paired end?
      pe=1
      shift 1
      ;;
    -d | --dir) #directory of inputs
      indir=$2
      shift 2
      ;;
    -s | --samp) #sample basename
      samp=$2
      shift 2
      ;;
    *) break
      ;;
  esac
done
## inputs
read1s=`ls ${indir}/*${samp}*1.fastq.gz`

###trimming and alignment
if [ -z $pe ]
then
  ## sing-end reads
  echo "single-end input"
  echo "trimming $read1s"
  parallel --no-notice ~/Code/trim_galore_zip/trim_galore -q 28 \
   -o ${outdir} {} ::: $read1s
#  ~/Code/trim_galore_zip/trim_galore -q 28 \
#    ${reads} \
#    -o ${outdir} \
#    &> ${outdir}/${samp}_trimlog
else
  ## paired-end reads
  echo "paired-end input"
  echo "trimming $read1s"
  ## code for getting the read pairs from read1s
  parallel --no-notice ~/Code/trim_galore_zip/trim_galore -q 28 \
    -o ${outdir} \
    --paired {} {=s/_1.fastq_/_2.fastq_/=} ::: ${read1s}
fi
