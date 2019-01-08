#!/bin/bash
dmdir=/shared/Data/demux
[ -e $dmdir ]||mkdir $dmdir
indir=/shared/raw

fqs=$(find $indir -maxdepth 1 -name "*f*q.gz")
for fq in $fqs; do
  base=$(basename "$fq")
  base=${base%%.*}
  outdir=$dmdir/$base
  [ -e $outdir ]||mkdir -p $outdir
#  bases="$bases $base"
  log=$indir/$base.demux.log
  com="porechop --threads 72 -i $fq -b $outdir &> $log"
  echo $com
  eval $com
done
#parallel -j 1 $com ::: $bases
