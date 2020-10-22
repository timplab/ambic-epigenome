#!/bin/bash
root=/shared/Data
bamdir=$root/bam

samples=$(awk 'NR>1' $root/chosigma_samples.txt | tr "\n" " ")
for samp in $samples; do
  outbam=$bamdir/$samp.chok1gshd_sigmaIgG.bam
  bams=$(find $bamdir -name "*$samp*bam")
  if [[ $(echo $bams | wc -w) -gt 1 ]];then
    echo $bams
    samtools merge -f $outbam $bams
    samtools index $outbam
  else
    bai=$(readlink -f $bamdir/$samp*.bai)
#    echo $bai
#    mv $bams $outbam
#    mv $bai "$outbam.bai"
  fi
done
#parallel samtools index $bamdir/{}.chok1gshd_sigmaIgG.bam ::: $samples
