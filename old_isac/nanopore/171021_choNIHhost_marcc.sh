#!/bin/bash
root=/scratch/users/ilee29@jhu.edu
tarpath=$root/tar
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
#listpath=$root/choNIHhost_list.txt

for raw in `ls $tarpath/*cho*`;do
  base=$(basename "$raw")
  base=${base%%.raw*}
  echo $base
  wpath=$root/$base
  mkdir -p $wpath
  tarpath=$root/tar
  tar=`ls ${tarpath}/${base}.raw.tgz` 
  rawpath=${wpath}/raw
  callpath=${wpath}/called
  extpath=${wpath}/ext
  fast5path=${wpath}/called_fast5
  mkdir $rawpath
  
#  echo "untarring $tar"
#  sbatch -o ../untar_$base.log -e ../untar_$base.log -D $rawpath \
#    $slurmpath/untar.scr $tar 
  
#  # identify fast5 folder path and perform basecalling
#  fast5path=`readlink -f ${rawpath}`
#  echo "basecalling reads in $fast5path"
#  $slurmpath/call_wrapper.sh -i $fast5path -w $wpath \
#    -f FLO-MIN106 -s $slurmpath --marcc &> $wpath/callwrapper.log
#  #
  
  echo "concatentating the fastq files"
  sbatch -D $wpath --partition=shared $slurmpath/cat_bcall.scr -w $wpath -p ${base}
  
done
