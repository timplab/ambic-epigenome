#!/bin/bash
root=/scratch/users/ilee29@jhu.edu/171025_cho
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
logpath=$root/log
rawtarpath=$root/rawtar

for raw in `ls $rawtarpath/* `;do
  base=$(basename "$raw")
  base=${base%%.*}
  echo $base
  rawpath=$root/raw/$base
  if [ -e $rawpath ];then
    echo "already done, moving on"
  else
    echo "extracting $raw to $rawpath"
    mkdir $rawpath
    sbatch -o $logpath/untar_$base.log -e $logpath/untar_$base.log -D $rawpath \
      $slurmpath/untar.scr $raw 
  fi
  
  
#  # identify fast5 folder path and perform basecalling
#  fast5path=`readlink -f ${rawpath}`
#  echo "basecalling reads in $fast5path"
#  $slurmpath/call_wrapper.sh -i $fast5path -w $wpath \
#    -f FLO-MIN106 -s $slurmpath --marcc &> $wpath/callwrapper.log
#  #
  
#  echo "concatentating the fastq files"
#  sbatch -D $wpath --partition=shared $slurmpath/cat_bcall.scr -w $wpath -p ${base}
  
done
