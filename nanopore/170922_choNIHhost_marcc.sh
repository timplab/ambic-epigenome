#!/bin/bash
root=/scratch/users/ilee29@jhu.edu
srcpath=/home-2/ilee29@jhu.edu/Code/ilee/oxford
slurmpath=$srcpath/slurm
base1=170922_choNIHhost1
base2=170922_choNIHhost2

for base in $base1 $base2;do
wpath=$root/$base
mkdir -p $wpath
tarpath=$wpath/tar
tar=`ls ${tarpath}/${base}.raw.tgz` 
rawpath=${wpath}/raw
callpath=${wpath}/called
extpath=${wpath}/ext
fast5path=${wpath}/called_fast5

#echo "untarring $tar"
#sbatch -o untar_$base.log -e untar_$base.log -D $wpath \
#  $slurmpath/untar.scr $tar 

## identify fast5 folder path and perform basecalling
#fast5path=${rawpath}
#echo "basecalling reads in $fast5path"
#$slurmpath/call_wrapper.sh -i $fast5path -w $wpath \
#  -f FLO-MIN107 -s $slurmpath --marcc &> $wpath/callwrapper.log
##

#echo "extracting the fastq sequences"
#$slurmpath/extract_wrapper.sh -s $slurmpath -w $wpath --marcc > $wpath/extwrapper.log

echo "concatentating the fastq files"
sbatch -D $wpath --partition=shared $slurmpath/cat_extract.scr -w $wpath -p ${base}

done

#echo "archiving called reads in $callpath"
#$slurmpath/splitArchive.sh -i $callpath -n 5 -p ${tarpath}/${base}.called


