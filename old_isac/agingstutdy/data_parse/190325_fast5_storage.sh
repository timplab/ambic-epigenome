#!/bin/bash
s3root="s3://timp.ambic/choSigmaAgingStudy/nanopore/reads/"
for dir in $(aws s3 ls $s3root ); do
  if [ "$dir" == "PRE" ];then
    continue
  fi
  run=${dir%/}
  size=$(aws s3 ls $s3root$dir | awk '{ sum+=$3}END{ print sum }')
  echo $run,$(($size/1000000000))
done
