#!/bin/bash
TMPDIR=/mnt/c/tmp
export TMPDIR=$TMPDIR TMP=$TMPDIR TEMP=$TMPDIR

root=/mnt/c
outdir=$root/pooled/fastq
days="0 30 60 90"
for day in $days; do
  for cell in Host StableGln StableNogln UnstableGln UnstableNogln; do
    for rep in 1 2 3; do
      base="CHOZN${cell}Day${day}_$rep"
      bases="$bases $base"
    done
  done
done
com="find /mnt/*/reads -maxdepth 1 -type f -name \"{}*fastq.gz\" \
  -exec gunzip -c \{\} \\; | gzip > $outdir/{}.fastq.gz"
#  com="echo {}"
#parallel -j 15 $com ::: $bases

# sanity check
#  com="gunzip -c $outdir/{}.fastq.gz | awk 'NR%4==1{ print }' > $outdir/{}.fqnames.txt"
#  com="grep sampleid $outdir/{}.fqnames.txt | wc -l > $outdir/{}.fqnum.txt"

parallel $com ::: $bases
