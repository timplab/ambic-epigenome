#!/bin/bash
root=/dilithium/Data/NGS/projects/choSigmaAgingStudy
# change these
fqdir=$root/run2/fastq
[ -e $fqdir ]||mkdir $fqdir
runid=wtimp1_181126_choatac_2

if [ "$1" == "dl" ];then
  ascp -v cidr-timp@162.129.245.24:$runid/ $root
fi

if [ "$1" == "organize" ];then
  csv=$root/$runid/$runid.csv
  indir=$root/$runid/FASTQ
  awk 'NR>1' $csv | while IFS=$',' read -r -a line
  do
    if [ "${line[2]}" -eq 2 ];then
      continue
    fi
    idx=${line[3]}
    name=${line[4]}
    echo $name
    out=$fqdir/$name.fastq.gz
    fqs=$(readlink -f "$indir/*$idx*fastq.gz")
    zcat $fqs | gzip > $out
  done
fi

if [ "$1" == "merge" ];then
  sampsfh=$root/181129_choagingstudy_atacseq_samples.txt
  fqdir=$root/fastq
  [ -e $fqdir ]||mkdir $fqdir
  samps=$(awk 'NR>1{print}' $sampsfh | tr '\r\n' ' ')
  for samp in $samps; do
    echo $samp
    outfq=$fqdir/$samp.fastq.gz
    ls $root/run*/fastq/*$samp.fastq.gz
    cat $root/run*/fastq/*$samp.fastq.gz > $outfq
  done
fi
