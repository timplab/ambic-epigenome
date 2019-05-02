#!/bin/bash

## $1 = Directory containing all vcf files
## $2 = Output directory

ls ${1}*.vcf > ${1}vcf_files.txt

echo Merging files: `cat ${1}/vcf_files.txt`

SURVIVOR merge ${1}/vcf_files.txt 1000 1 1 -1 -1 -1 $2/merged_SURVIVOR_1kbpdist_typesave.vcf
