#!/bin/bash
root=/kyber/Data/Nanopore/projects/ambic/sigma
bdir=$root/bed
cont=$bdir/sigma_in_host.bed
rnames=$bdir/sigma_in_host_rnames.txt
runs=$bdir/sigma_in_host_runs.txt
dbdir=$root/readdb

if [ "$1" == "find" ];then
#  grep sigma $bdir/* > $cont
  echo $(cut -f4 $cont | sort | uniq -c | awk '{ print $2 }') | tr " " "\n" > $rnames
fi

if [ "$1" == "run" ];then
  grep -f $rnames $dbdir/* > $runs
fi




