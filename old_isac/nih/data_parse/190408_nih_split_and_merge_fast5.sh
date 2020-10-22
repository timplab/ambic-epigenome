#!/bin/bash
rdir=/mnt/c/Data/reads
mergedir=/mnt/e/Data/reads_merged
t=72
for run in $(find $rdir/* -maxdepth 0 -type d); do
  samp=$(basename $run)
  s3num=$(aws s3 ls s3://timp.ambic/nih/nanopore/reads/$samp/ | grep fast5 | wc -l)
  if [ $s3num -gt 0 ];then
    echo $samp,$s3num
    continue
  else
    echo $run
  fi
  if [ "$1" == "dryrun" ];then
    continue
  fi
  outdir=$mergedir/$samp
  [ ! -e $outdir ]||rm -r $outdir
  mkdir -p $outdir
  file1=$(find $run/* -maxdepth 0 -print -quit)
  if [[ $file1 == *fast5 ]]; then
    # first split 
    echo "splitting"
    dir=$dev/reads_intermediate/$samp
    [ ! -e $dir ]||rm -r $dir
    mkdir -p $dir
    multi_to_single_fast5 --recursive -t $t \
      -i $run -s $dir
  else
    dir=$run
  fi
  echo "merging"
  log=$outdir.merge.log
  single_to_multi_fast5 --recursive -t $t -n 4000 \
    -f $samp -i $dir -s $outdir &> $log
  echo "syncing"
  aws s3 sync $outdir s3://timp.ambic/nih/nanopore/reads/$samp
  echo "removing"
  rm -r $outdir
  if [[ $file1 == *fast5 ]]; then
    rm -r $dir
  fi
done
