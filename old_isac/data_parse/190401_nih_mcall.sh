#!/bin/bash
root=/mnt/c/Data
TMP=$root/tmp
TEMP=$TMP
TMPDIR=$TMP
export TMP TEMP TMPDIR
mdir=$root/mcall
rdir=$root/reads
bdir=$root/bam
mbdir=$root/mbed
[ -e $mbdir ]||mkdir $mbdir
mfdir=$root/mfreq
[ -e $mfdir ]||mkdir $mfdir
ref=$root/Reference/Cricetulus_griseus_chok1gshd.nihIgG.fa
np=/home/ubuntu/Code/nanopolish/nanopolish

if [ "$1" == "index" ];then
  for dir in $(find $rdir/* -maxdepth 0 -type d); do
    samp=$(basename "$dir")
    samps="$samps $samp"
  done
  dir=$rdir/{}
  log=$dir.index.log
  fq=$dir.fastq.gz
  sum=$dir.summary.txt
  com="$np index --verbose -s $sum -d $dir $fq &> $log"
  parallel $com ::: $samps
fi

if [ "$1" == "mcall" ];then
  for dir in $(find $rdir/* -maxdepth 0 -type d); do
    samp=$(basename "$dir")
    fq=$dir.fastq.gz
    bam=$bdir/$samp.sorted.bam
    out=$mdir/$samp.cpg.meth.tsv.gz
    log=$mdir/$samp.cpg.mcall.log
    $np call-methylation --verbose -g $ref -r $fq -b $bam -t 72 2> $log |\
      pigz -p 18 -c > $out
  done
fi

if [ "$1" == "mbed" ];then
  scr=/home/ubuntu/Code/ilee/nanopolish/mtsv2bedGraph.py
  for tsv in $(find $mdir -name "*tsv.gz"); do
    base=$(basename "$tsv")
    base=${base%%.*}
    bases="$bases $base"
    tmpdir=$root/tmp/$base
    [ -e $tmpdir ]||mkdir $tmpdir
  done
  mtsv="$mdir/{}.cpg.meth.tsv.gz"
  mbed="$mbdir/{}.cpg.meth.bed.gz"
  com="gunzip -c $mtsv | python $scr |\
    sort -T $root/tmp/{} -k1,1 -k2,2n | bgzip > $mbed \
    && tabix -p bed $mbed"
#  echo $com
#  com="echo {}"
  parallel $com ::: $bases
fi

if [ "$1" == "mfreq" ];then
  scr=/home/ubuntu/Code/nanoNOMe/scripts/parseMethylbed.py
  for tsv in $(find $mdir -name "*tsv.gz"); do
    base=$(basename "$tsv")
    base=${base%%.*}
    bases="$bases $base"
    tmpdir=$root/tmp/$base
    [ -e $tmpdir ]||mkdir $tmpdir
  done
  samp="{}"
  mbed="$mbdir/{}.cpg.meth.bed.gz"
  out="$mfdir/{}.cpg.meth.freq.txt.gz"
  com="python $scr frequency -i $mbed |\
    bgzip > $out && tabix -b 2 -e 2 $out"
#  echo $com
  parallel $com ::: $bases
fi


