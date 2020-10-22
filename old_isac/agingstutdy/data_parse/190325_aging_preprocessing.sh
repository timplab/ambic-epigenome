#!/bin/bash
root=/mnt
f5root=$root/c/reads
datroot=$root/d/Data
fqdir=$datroot/fastq
bamdir=$datroot/bam
mdir=$datroot/mcall
refpre=$datroot/ref/Cricetulus_griseus_chok1gshd.sigmaIgG
mmi=$refpre.mmi
fa=$refpre.fa
np=~/Code/nanopolish/nanopolish
awsroot="s3://timp.ambic/choSigmaAgingStudy/nanopore"
plog=$datroot/process.log
echo "sample,fqsize(kb),time(s)" > $plog

for fq in $(aws s3 ls $awsroot/fastq/ | grep -v "index" | grep Day0 | awk '{ print $NF }'); do
  run=${fq%%.*}
  f5dir=$f5root/$run
  start=$(date +%s)
  if [ $run == "CHOZNHostDay30_2_PAD05424" ];then 
    continue
  fi
  if [ "$1" == "clean" ];then
    rm -r $f5root/* $fqdir/* $bamdir/* $mdir/*
    exit
  fi
  maws=$(aws s3 ls $awsroot/mcall/ | grep $run)
  if [ -n "$maws" ];then
    echo $run already done - skipping
    continue
  fi
  echo $run
  if [ "$1" == "dryrun" ];then
    continue
  fi
  # first download
  [ -e $f5dir ]||mkdir $f5dir
  bam=$bamdir/$run.sorted.bam
  log=$bamdir/$run.align.log
  indexlog=$fqdir/$run.index.log
  # align and dl f5/index in parallel
  aws s3 sync $awsroot/fastq/ $fqdir --exclude "*" --include "$run*" --quiet &
  aws s3 sync $awsroot/bam/ $bamdir --exclude "*" --include "$run*" --quiet &
  aws s3 sync $awsroot/reads/$run/ $f5dir --quiet &
  wait 
  echo "indexing"
  $np index --verbose -d $f5dir $fqdir/$fq &> $indexlog
  fqsize=$(du $fqdir/$fq | cut -f1)
  # call methylation
  mout=$mdir/$run.cpg.meth.tsv.gz
  echo calling methylation
  log=$mdir/$run.mcall.log
  $np call-methylation --verbose -t 72 -r $fqdir/$fq -b $bam -g $fa 2> $log |\
    gzip > $mout
  # upload everything
  echo uploading
  aws s3 sync $fqdir $awsroot/fastq/ --quiet
  aws s3 sync $bamdir $awsroot/bam/ --quiet
  aws s3 sync $mdir $awsroot/mcall/ --quiet
  # delete
  rm $f5dir/* $fqdir/* $bamdir/* $mdir/*
  end=$(date +%s)
  elapse=$((${end}-${start}))
  echo $run,$fqsize,$elapse >> $plog
done
