#!/bin/bash
root=/dilithium/Data/Nanopore/Analysis/171024_cho
rawroot=/dilithium/Data/Nanopore/oxford
base=choNIHhost

for i in 1 2 3;
do
  rep=${base}$i
  echo $rep
  find $rawroot/*$base* -name "*$rep*fq.gz" -exec cat {} \; > $root/fastq/$rep.fq.gz
done

