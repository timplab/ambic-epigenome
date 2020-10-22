#!/bin/bash
root=/scratch/users/ilee29@jhu.edu
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
listpath=$root/choNIHhost_list.txt

for base in `cat $listpath`;do
  echo $base
  wpath=$root/$base
#  mkdir -p $wpath
  tarpath=$root/tar
#  tar=`ls ${tarpath}/${base}.raw.tgz` 
  rawpath=${wpath}/raw
  callpath=${wpath}/called
  extpath=${wpath}/ext
  fast5path=${wpath}/called_fast5
  
#  echo "untarring $tar"
#  sbatch -o untar_$base.log -e untar_$base.log -D $wpath \
#    $slurmpath/untar.scr $tar 
  
#  # identify fast5 folder path and perform basecalling
#  fast5path=`readlink -f ${rawpath}`
#  echo "basecalling reads in $fast5path"
#  $slurmpath/call_wrapper.sh -i $fast5path -w $wpath \
#    -f FLO-MIN106 -s $slurmpath --marcc &> $wpath/callwrapper.log
#  #
  
#  echo "extracting the fastq sequences"
#  $slurmpath/extract_wrapper.sh -s $slurmpath -w $wpath --marcc > $wpath/extwrapper.log
  
  echo "concatentating the fastq files"
  sbatch -D $wpath --partition=shared $slurmpath/cat_extract.scr -w $wpath -p ${base}
  
done

#echo "archiving called reads in $callpath"
#$slurmpath/splitArchive.sh -i $callpath -n 5 -p ${tarpath}/${base}.called


