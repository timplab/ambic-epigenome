#!/bin/bash
svout=/dilithium/Data/Nanopore/Analysis/171025_cho/sniffles/choNIHIgG.sniffles.vcf
grep 4_0 $svout | awk 'FS=";"{ print $3 }' | grep -v 4_0 | uniq
grep 4_0 $svout | cut -f1,2,3,10 | grep -v 4_0 | uniq
