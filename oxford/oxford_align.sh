#!/bin/bash

while :
do
  case "$1" in
    -r | --reference) #reference
      ref=$2
      shift 2
      ;;
    -i | --input) #input
      in=$2
      shift 2
      ;;
    -a | --aligner) #which aligner?
      aligner=$2
      shift 2
      ;;
    -o | --outdir) #output dir; default: input dir
      outdir=$2
      shift 2
      ;;
    -s | --sample) #sample name to be used as prefix
      samp=$2
      shift 2
      ;;
    --samarg) #extra argument for samtools
      samarg=$2
      shift 2
      ;;
    *) break
      ;;
  esac
done

srcdir=/home/isac/Code
if [ $outdir ];then
  outpre=${outdir}/${samp}
else
  outpre=$(dirname "$input")/${samp}
fi
echo "output prefix: $outpre"

case "$aligner" in
  ngmlr)
    echo "using ngmlr to align"
    args="-t 10 -r ${ref} -q ${in} -x ont"
    src=${srcdir}/ngmlr/bin/ngmlr-0.2.6/ngmlr
    ;;
  bwa)
    echo "using bwa mem to align"
    args="-t 10 -x ont2d ${ref} ${in}"
    src="bwa mem"
    ;;
  minimap2)
    echo "using minimap2 to align"
    refind=${ref%.*}.mmi
    if [ ! -e $refind ];then
      echo "make index first with : minimap2 -d ref.mmi ref.fa"
      exit
    fi
    args="-ax map-ont $refind $in"
    src=${srcdir}/minimap2/minimap2
    ;;
esac


echo "$src $args"
$src $args | \
  samtools view $samarg -Sb - | \
  samtools sort -o ${outpre}.sorted.bam
samtools index ${outpre}.sorted.bam

## index stats
samtools idxstats ${outpre}.sorted.bam > ${outpre}_idxstats.txt
awk '{ sum+= $3 }END { print "aligned =",sum,"unaligned =",$4,"fraction =",sum/(sum+$4) }' ${outpre}_idxstats.txt
