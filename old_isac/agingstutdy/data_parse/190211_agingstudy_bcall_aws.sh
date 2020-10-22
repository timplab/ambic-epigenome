#!/bin/bash
root=/shared
readroot=$root/reads
bcallroot=$root/bcall
[ -e $bcallroot ]||mkdir $bcallroot
outdir=/shared/fastq
np=/home/ubuntu/Code/nanopolish/nanopolish

for dir in $(find $readroot/* -maxdepth 0 -type d); do
	samp=$(basename "$dir")
  samps="$samps $samp"
done

if [ "$1" == "bcall" ];then
  src=/home/ubuntu/Code/ilee/oxford/slurm/call_wrapper.sh
  for samp in $samps; do
    cfg=~/Code/ilee/oxford/config/guppy_r941_10kchunk.cfg
    if [[ $samp == *"PAD"* ]];then
      cfg=~/Code/ilee/oxford/config/guppy_r941_10kchunk_prom.cfg
    fi
    indir=$readroot/$samp
    outdir=$bcallroot/$samp
    $src -i $indir -o $outdir -a "-c $cfg" -e aws --guppy 
  done
fi

if [ "$1" == "cat" ];then
  indir="$bcallroot/{}"
  fqs="$indir/*/*fastq"
  sums="$indir/*/sequencing_summary.txt"
  outfq=$outdir/{}.fastq.gz
  outsum=$outdir/{}.summary.txt
  echo sum
  com="awk 'NR==1||\$1!=\"filename\"{ print }' $sums > $outsum"
  parallel $com ::: $samps
  echo fq
  com="cat $fqs | gzip > $outfq"
  parallel $com ::: $samps
fi

if [ "$1" == "check" ];then
  fq=$outdir/{}.fastq.gz
  sum=$outdir/{}.summary.txt
  out=$outdir/{}.fqnames.txt
  com="gunzip -c $fq | awk 'NR%4==1{print}' > $out"
#  parallel $com ::: $samps
  nums=$outdir/fqnum.txt
#  [ -e $nums ]&&rm $nums
#  for samp in $samps; do
#    out=$outdir/$samp.fqnames.txt
#    num=$(grep sampleid $out | wc -l)
#    echo $num $samp >> $nums
#  done
  sumnum=$outdir/summarynum.txt
#  wc -l $outdir/*summary.txt > $sumnum
  ./compare_readnum.py $nums $sumnum
fi

if [ "$1" == "index" ];then
  outdir=$readroot
  fq=$outdir/{}.fastq.gz
  dir=$outdir/{}
  sum=$outdir/{}.summary.txt
  log=$outdir/{}.index.log
  com="$np index -v -d $dir -s $sum $fq &> $log"
  parallel --jobs 72 $com ::: $samps
fi
