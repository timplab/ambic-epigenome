#!/bin/bash
root=/kyber/Data/Nanopore/projects/ambic/capture/190225_choSigma_plasmid_target_pipeline_test/data/nanopore
rdir=$root/reads
bamdir=$root/bam
mdir=$root/mcall
bdir=$root/mbed
idir=$root/igv
ref=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.fa

cd $rdir

for samp in $(find * -maxdepth 0 -type d); do
  fq=$samp.fastq.gz
  db=$fq.index.readdb
  if [ ! -e $db ];then
    echo $samp
    nanopolish index -d $samp $fq
  fi
  bam=$bamdir/$samp.minimap2.sorted.bam
  if [ ! -e $bam ];then
    echo "aligning $samp"
    log=$bamdir/$samp.minimap2.align.log
    idx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.mmi
    minimap2 --MD -t 10 -ax map-ont -L $idx $fq 2> $log |\
      samtools view -b -q 20 - |\
      samtools sort -o $bam
    samtools index $bam
  fi
  mtsv=$mdir/$samp.cpg.meth.tsv.gz
  if [ ! -e $mtsv ];then
    echo "calling methylation on $samp"
    log=$mdir/$samp.mcall.log
    nanopolish call-methylation --verbose -t 10 \
      -r $fq -b $bam -g $refo 2> $log |\
      gzip > $mtsv
  fi
  # make mbed
  mbed=$bdir/$samp.cpg.meth.bed.gz
  if [ ! -e $mbed ];then
    echo "making $samp mbed"
    scr=/home/isac/Code/nanopore-methylation-utilities/mtsv2bedGraph.py
    log=$bdir/$samp.mbed.log
    $scr -i $mtsv | sort -k1,1 -k2,2n | bgzip > $mbed
    tabix -p bed $mbed
  fi
  # igv
  mbam=$idir/$samp.methylation.bam
  if [ "$1" == "mbam" ];then
    echo "converting $samp bam"
    scr=/home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py
    log=$idir/$samp.convert.log
    $scr -t 10 --verbose -c $mbed -b $bam 2> $log |\
      samtools sort -o $mbam
    samtools index $mbam
  fi
done

