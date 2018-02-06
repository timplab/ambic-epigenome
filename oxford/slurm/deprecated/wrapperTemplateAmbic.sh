#!/bin/bash
root=/shared/data
base=170710_corden_WT
base=170825_corden_mut
echo $base
#tar=`ls ${tarpath}/${base}.raw.tgz` 
wpath=/shared/data/$base
mkdir -p $wpath
tarpath=$wpath/tar
srcpath=/shared/Code/ilee/oxford/slurm
slurmpath=/shared/Code/ilee/slurm
rawpath=${wpath}/raw
callpath=${wpath}/called
extpath=${wpath}/ext
fast5path=${wpath}/called_fast5

#echo "untarring $tar"
#sbatch -o untar_$base.log -e untar_$base.log -D $wpath \
#  $srcpath/untar.scr $tar 

## identify fast5 folder path and perform basecalling
#fast5path=${rawpath}
#echo "basecalling reads in $fast5path"
#$srcpath/call_wrapper.sh -i $fast5path -w $wpath \
#  -f FLO-MIN107 -s $srcpath 
##


#echo "extracting the fastq sequences"
#$srcpath/extract_wrapper.sh -s $srcpath -w $wpath 


#echo "concatentating the fastq files"
#sbatch -D $wpath $srcpath/cat_extract.scr -w $wpath -p ${base}

#
echo "archiving called reads in $callpath"
$slurmpath/splitArchive.sh -i $callpath -n 5 -p ${tarpath}/${base}.called


