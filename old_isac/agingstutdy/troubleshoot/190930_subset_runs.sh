#!/bin/bash
dir="/kyber/Data/Nanopore/projects/ambic/sigma/methylation/mbed"
outdir=$dir/../subset_troubleshoot
[ -e $outdir ]||mkdir $outdir

if [ "$1" == "subset" ];then
   du -s $dir/*bed.gz  | awk '$1>000000{print $2}' | xargs -d "\n" mv -t $outdir
fi

if [ "$1" == "mfreq" ];then
  for mtsv in $(find $outdir -name "*meth.bed.gz"); do
    base=$(basename "${mtsv%%.bed.gz}")
    bases="$bases $base"
  done
  base="{}"
  mtsv=$outdir/$base.bed.gz
  out=$outdir/$base.freq.txt.gz
  scr="/home/isac/Code/nanopore-methylation-utilities/parseMethylbed.py frequency"
  com="gunzip -c $mtsv | $scr | bgzip > $out"
  parallel $com ::: $bases
fi

