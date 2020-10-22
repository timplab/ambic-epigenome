#!/bin/bash
rawdir=/kyber/Data/NGS/Raw/190310_chosigmaday0_atacseq/wtimp1_154135/FASTQ
root=/kyber/Data/NGS/projects/ambic/agingstudy/atacseq/190310_chosigmaday0_atacseq/data/atacseq
fqdir=$root/fastq
codedir=/home/isac/Code/projects/ambic-epigenome
if [ "$1" == "org" ];then
  csv="/home/isac/Dropbox/Data/ambic/aging_study/atac-seq/190313_atacseq_day0_samples.csv"
  while IFS=$',' read -r -a line || [[ -n "$line" ]]
  do
    samp=${line[0]}
    bc=$(echo ${line[1]} | tr -d "\r")
    fqs=$(find $rawdir -name "*$bc*fastq.gz")
    args="$args $bc $samp"
  done <<< "$(awk 'NR>1' $csv)"
  com="zcat $rawdir/*{1}*fastq.gz | gzip > $fqdir/{2}.fastq.gz"
#  com="echo {1} {2}"
  parallel -N 2 $com ::: $args
fi

cd $codedir
snakemake --cores 10 parse_atacseq
