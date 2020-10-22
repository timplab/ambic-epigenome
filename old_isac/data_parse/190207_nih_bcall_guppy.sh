#!/bin/bash
root=/shared/Data
f5root=$root/reads
bcallroot=$root/bcall
outdir=$f5root
cd $f5root

for samp in $(find * -maxdepth 0 -type d); do
  samps="$samps $samp"
done
if [ "$1" == "bcall" ];then
  cfg=~/Code/ilee/oxford/config/guppy_10kchunk.cfg
  for samp in $samps; do
    f5dir=$f5root/$samp
    bcalldir=$bcallroot/$samp
    [ -e $bcalldir ]||mkdir $bcalldir
    for dir in $(find $f5dir/* -maxdepth 0 -type d); do
      ind=$(basename "$dir")
      outdir=$bcalldir/$ind
      log=$bcalldir/$ind.bcall.log
      if [ -e $log ];then
        msg=$(tail -n1 $log)
        if [ "$msg" == "Basecalling completed successfully." ];then
          continue
        else
          rm $log
          rm -r $outdir
        fi
      fi
      echo $outdir
      com="guppy_basecaller --verbose_logs -x cuda:0 \
        -i $dir --num_callers 16  -s $outdir --recursive -q 0 -c $cfg &>> $log"
      echo $com > $log
      eval $com
    done
  done
fi

if [ "$1" == "cat" ];then
  bcalldir=$bcallroot/{}
  fqs="$bcalldir/*/*fastq"
  outfq="$outdir/{}.fastq.gz"
  echo "fq"
  com="cat $fqs | gzip > $outfq"
  parallel $com ::: $samps
  echo "sum"
  for samp in $samps; do
    sums="$bcallroot/$samp/*/sequencing_summary.txt"
    outsum="$outdir/$samp.summary.txt"
    awk 'NR==1||$1!="filename"{ print}' $sums > $outsum
  done
fi

if [ "$1" == "getnum" ];then
  for cell in host IgG; do
    for rep in 1 2 3; do
      sums=$(find $outdir -maxdepth 1 -name "*$cell$rep*summary.txt")
      totbp=$(awk '{sum+=$12}END{print sum}' $sums)
      echo $cell$rep,$totbp
    done
  done
fi

if [ "$1" == "merge" ];then
  for cell in host IgG; do
    for rep in 1 2 3; do
      base=choNIH$cell$rep
      bases="$bases $base"
      sums=$(find $f5root/$base -maxdepth 1 -name "$base*summary.txt")
      outsum=$f5root/$base.summary.txt
      awk 'NR==1||$1!="filename"{ print}' $sums > $outsum
    done
  done
  com="zcat $f5root/{}/{}*fastq.gz | gzip > $f5root/{}.fastq.gz"
  parallel $com ::: $bases
fi

