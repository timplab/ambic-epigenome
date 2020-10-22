#!/bin/bash
indir=/data/basecalled
outroot=/data/demux
samps="190314_choSigmaHostD30r3StaGlnD90r2NoglnD30r3 190314_choSigmaUnstaGlnD60r2D90r1NoglnD0r3"
for samp in $samps; do
  echo $samp
  if [ "$1" == "demux" ];then
    outdir=$outroot/$samp
    log=$outroot/$samp.demux.log
    guppy_barcoder --recursive -t 90 -i $indir/$samp -s $outdir > $log
  fi
done

