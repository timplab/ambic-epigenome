#!/bin/bash

## First sniffles command. Takes single .bam file as input
# Directory has to be structured: 
#	$WD/log/
#	$WD/vcf/1-originalVcf/
#	$WD/tmp/

## Reads supporting SV = 2 (-s 2)
ml gcc/5.5.0

WD=/home-3/kmcfarl6@jhu.edu/scratch/epigenomicsAMBIC/nanopore/svAnalysis/sigma
FILE=$1
BASE=`basename $FILE`
LOGFILE=${WD}/log/rawCalled.${BASE}
OUTFILE=$WD/vcf/1-originalVcf/${BASE::-4}.vcf
TEMPFILE=$WD/tmp/rawCalled.${BASE}

echo "Creating $OUTFILE from $FILE"
sniffles -t $(nproc) -s 2 -n -1 --cluster --tmp_file ${TEMPFILE}.tmp -m ${1} -v ${OUTFILE}.vcf &> ${LOGFILE}.log
echo "Finished with $OUTFILE"


