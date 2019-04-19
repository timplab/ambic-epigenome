#!/bin/bash

## $1 = Directory containing all vcf files
## $2 = Output directory

for file in ${1}*.vcf; do
	echo Sorting ${file}
	bcftools sort -Ov -o ${2}/`basename $file`.sorted $file
	echo $file sorted
done

rename .vcf.sorted .sorted.vcf ${2}/*
