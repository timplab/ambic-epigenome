#!/bin/bash

## Second sniffles command. Force call command script. $1 is individual .bam file being evaluated
## Directory has to be structured:
#	$WD/log/
#	$WD/vcf/4-forceCalled/
#	$WD/tmp/
#
# Reads supporting SV = 5 (-s 5) 

ml gcc/5.5.0

WD=/home-3/kmcfarl6@jhu.edu/scratch/epigenomicsAMBIC/nanopore/svAnalysis/sigma
FILE=$1
BASE=`basename $FILE`
LOGFILE=${WD}/log/forceCalled.${BASE}
OUTFILE=${WD}/vcf/4-forceCalled/${BASE::-4}
TEMPFILE=${WD}/tmp/forceCalled.${BASE}
MERGEDLIST=${WD}/vcf/3-mergedOrig/merged_SURVIVOR_1kbpdist_typesave.vcf


echo "Force-calling $OUTFILE with $1 using ${MERGEDLIST}"
sniffles -t $(nproc) -s 5 -n -1 --Ivcf ${MERGEDLIST} --cluster --tmp_file ${TEMPFILE}.tmp -m ${1} -v ${OUTFILE}.fc.vcf &> ${LOGFILE}.log && echo "Done with ${OUTFILE}.vcf"
