#!/bin/bash
srcdir=$(dirname "$0")
reorgscr=${srcdir}/../../../oxford/combinatar.py

dirs="/mnt/ambic3/reads /mnt/ambic-fast5/reads"

for dir in $dirs; do
  runs=$(find $dir/unorganized/* -maxdepth 0 -type d)
  for run in $runs; do
    base=$(basename "$run")
    echo $base
    python -u $reorgscr -n 100000 -s $base -o $dir -i $run #--dryrun
  done
done



