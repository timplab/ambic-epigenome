#!/bin/bash

## Usage: ./atacAlign.sh -o [outdir] -x [index path] \
##  -s [optional; sample name] [read1]|[read2]
# automatically determines PE or SE by number of input reads

###argument parsing
while :
do
  case "$1" in
    -i | --indir) ##input dir
      indir=$2
      shift 2
      ;;
    -o | --outdir) # directory for final output
      mkdir -p $2
      outdir=`readlink -f $2`
      shift 2
      ;;
    -x | --bowtieindex) # reference genome bowtie index
      btidx=$2
      shift 2
      ;;
    -s | --sample) # name of the output
      samp=$2
      shift 2
      ;;
    -t | --threads)
      t=$2
      shift 2
      ;;
    -p | --paired-end) #paired-end?
      pe=1
      shift 1
      ;;
    *) break
      ;;
  esac
done
## inputs
if [ -z $t ];then
  t=10
fi
##final output
sorted=${outdir}/${samp}.sorted.bam
echo outputting sorted alignment files to: $sorted
##tempdir
tmpdir=/atium/Data/tmp/$samp
if [ -d $tmpdir ];then 
  rm -r $tmpdir 
fi
mkdir -p $tmpdir

if [ -z $pe ]
then
  ## sing-end reads
  echo "single-end input"
  echo "$samp alignment"
  read1=`readlink -f ${indir}/${samp}_trimmed.fq.gz`
  bowtie2 -k 4 -p $t -t --local \
    -x $btidx -U $read1 \
    -S ${tmpdir}/${samp}.sam \
    2> ${outdir}/${samp}_bowtie2.log 
else
  ## paired-end reads
  echo "paired-end input"
  echo "$samp alignment"
  read1=`readlink -f ${indir}/${samp}*1.fq.gz`
  read2=`readlink -f ${indir}/${samp}*2.fq.gz`
  bowtie2 -k 4 -X 2000 -p $t -t --local \
    -x $btidx -1 $read1 -2 $read2 \
    -S ${tmpdir}/${samp}.sam
    2> ${outdir}/${samp}_bowtie2.log 
fi
## sam to bam and sort
samtools view -Su ${tmpdir}/${samp}.sam | \
  samtools sort -T ${tmpdir}${samp}.sorted -o ${sorted} -
## bam conversion, sort, and index
samtools index $sorted

## flagstat for QC
echo "$samp flagstat QC"
samtools flagstat $sorted &> ${outdir}/${samp}_flagstats.txt

## remove temp dir
rm -r $tmpdir
echo "$samp finished alignment"
