#!/bin/bash
root=/scratch/users/ilee29@jhu.edu/180222_choSigmaD0
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
logroot=$root/log
rawtarpath=$root/rawtar
rawroot=$root/raw
bcallroot=$root/bcall
config=/home-2/ilee29@jhu.edu/Code/ilee/oxford/config/albacore_config_sv_marcc.cfg

if [ "$1" == "untar" ];then
  logdir=$logroot/untar
  [ -e $logdir ]||mkdir -p $logdir
  for raw in `find $rawtarpath -name "*tgz"`;do
    base=$(basename "$raw")
    base=${base%%.*}
    echo $base
    rawpath=$root/raw/$base
    if [ -e $rawpath ];then
      echo "already done, moving on"
    else
      echo "extracting $raw to $rawpath"
      mkdir -p $rawpath
      log=$logdir/$base.untar.log
      sbatch -o $log -e $log -D $rawpath -p shared -t 6:0:0\
        $slurmpath/untar.scr $raw 
    fi
  done
fi

if [ "$1" == "bcall" ];then
  logdir=$logroot/bcall
  [ -e $logdir ]||mkdir $logdir
  for fast5path in `find $rawroot/* -maxdepth 0 -type d`;do
    base=$(basename "$fast5path")
    bcalldir=$bcallroot/$base
    [ -e $bcalldir ]||mkdir -p $bcalldir
    log=$logdir/$base.bcall.log
    echo "basecalling reads in $base"
  $slurmpath/call_wrapper.sh -i $fast5path -o $bcalldir \
    -s $slurmpath --marcc -a "--config $config" &> $log
  done
fi
  
#  echo "concatentating the fastq files"
#  sbatch -D $wpath --partition=shared $slurmpath/cat_bcall.scr -w $wpath -p ${base}
#  
#done
