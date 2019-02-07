#!/bin/bash

while :
do
  case "$1" in
    -d | --dir ) #dir
      dir="$2"
      shift 2
      ;;
    -b | --base ) #base
      base=$2
      shift 2
      ;;
    * ) break
      ;;
  esac
done
base=${base%%.*}
tsvs=$(find $dir -name "*$base*.tsv")

awk 'NR==1||FNR!=1' $tsvs
