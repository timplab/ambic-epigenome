#!/bin/bash


## this uses the assign_multimappers python script from encode atac-seq
##pipeline
srcdir=$(dirname "$0")
assign_multimappers=$srcdir/assign_multimappers.py

###argument parsing
while :
do
  case "$1" in
    -o | --outdir) # directory for final output
      mkdir -p $2
      outdir=`readlink -f $2`
      shift 2
      ;;
    -s | --sample) # name of the output
      samp=$2
      shift 2
      ;;
    -pe | --paired-end) # flag for paired-end alignments
      pe=1
      shift 1
      ;;
    *) break
      ;;
  esac
done
## input
sorted=`readlink -f ${outdir}/${samp}.sorted.bam`
echo "taking sorted alignment file $sorted"
if [ -z $samp ]
then
  samp=${sorted##*/}
  samp=${samp%%.*}
fi
echo "sample name : $samp"

##tempdir
tmpdir=/atium/Data/tmp/$samp
if [ -d $tmpdir ]; then
  rm -r $tmpdir
fi
mkdir -p $tmpdir
## intermediate outputs
qsorted=${tmpdir}/${samp}.sorted.bam
multimapped=${tmpdir}/${samp}.multi.sam
tmpfilt=${tmpdir}/${samp}.tmpfilt.bam
filt=${tmpdir}/${samp}.filt.bam
dupmark=${tmpdir}/${samp}.dupmark.bam

## final output
nodup=${outdir}/${samp}.nodup.bam

##src
picard=/home/isac/Code/picard.jar

##sambamba options
samop="-t 8 -p"
sortop="-t 8 -m 20G --tmpdir=$tmpdir -p"

if [ -z $pe ];then
  echo "sing-end bam processing"
  ## preprocess
  sambamba sort $sortop ${sorted} -n -o ${qsorted}
  samtools view -h $qsorted | $assign_multimappers -k 4 > $multimapped
  samtools view -F 1804 -Su $multimapped > $tmpfilt
  sambamba sort $sortop $tmpfilt -o $filt
else
  echo "paired-end bam processing"
  ## preprocess
  tmpview=${tmpdir}/${samp}.tmpview.bam
  fixmate=${tmpdir}/${samp}.fixmate.bam
  samtools view -F 524 -f 2 -u ${sorted} > $tmpview
  sambamba sort $sortop -n $tmpview -o ${qsorted}
  samtoolv view -h $qsorted | \
  $assign_multimappers -k 4 --paired-end > $multimapped
  samtools fixmate -r $multimapped $fixmate
  samtools view -F 1804 -f 2 -u $fixmate  > $tmpfilt
  sambamba sort $sortop $tmpfilt -o $filt
fi

## picard duplication marking
echo "marking duplicates using picard version"
echo "`java -jar $picard MarkDuplicates --version`"
joptions="-Xmx8G -Xms256M -XX:ParallelGCThreads=8 -Djava.io.tmpdir=${tmpdir}"
java $joptions -jar $picard MarkDuplicates \
  INPUT=$filt OUTPUT=$dupmark \
  METRICS_FILE=${outdir}/${samp}_duplicate_metric.txt \
  VALIDATION_STRINGENCY=LENIENT \
  ASSUME_SORTED=true REMOVE_DUPLICATES=false \
  &> ${outdir}/${samp}_dup.log

## remove duplicates
echo "removing duplicates"
if [ -z $pe ]; then
  samtools view -F 1804 -b $dupmark > $nodup
  sambamba index $samop $nodup
  sambamba flagstat $samop $nodup > ${outdir}/${samp}_nodup_flagstat.txt
  # complexity test
  # output:
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs
  # [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab]
  # PBC2=OnePair/TwoPair
  bedtools bamtobed -i $dupmark | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
  	grep -v 'chrM' | sort | uniq -c | \
  	awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} \
  ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
  END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
  printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n" \
  ,mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' \
    > ${outdir}/${samp}_pbcqc.txt
else
  samtools view -F 1804 -f 2 -b $dupmark > $nodup
  sambamba index $samop $nodup
  sambamba flagstat $samop $nodup > ${outdir}/${samp}_nodup_flagstat.txt
  ## complexity test
  ## output:
  # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs
  # [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab]
  # PBC2=OnePair/TwoPair
  sambamba sort $sortop -n $dupmark -o $dupmark.tmp.bam
	bedtools bamtobed -bedpe -i $dupmark.tmp.bam | \
		awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
		grep -v 'chrM' | sort | uniq -c | \
		awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} \
($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n"\
,mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' \
> ${outdir}/${samp}_pbcqc.txt
fi

## remove temp dir
rm -r $tmpdir
echo "finished with duplicate removal of ${samp}"
