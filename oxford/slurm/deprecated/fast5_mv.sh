#!/bin/bash
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford/marcc

if [ -z $1 ]; then
  echo "provide the path of the called fast5s"
  exit
else
  path=`readlink -f $1`
fi
mkdir called_fast5

for i in `ls $path`; do
  mkdir called_fast5/$i
  echo "moving called/${i}"
  mv called/${i}/workspace/*/*fast5 called_fast5/$i/
done

  
