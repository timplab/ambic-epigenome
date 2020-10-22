#!/bin/bash
bamdir=/atium/Data/Nanopore/Analysis/170608_cho/rawhit/
analysisdir=$bamdir
#bam=${bamdir}/170608_choIgGNIH_hits.sorted.bam
ngmlrbam=/atium/Data/Nanopore/Analysis/170608_cho/rawhit/hits.bwa.sorted.bam
bam=$ngmlrbam
bpnames=${analysisdir}/bpNames.txt
tnames=${analysisdir}/bptnames.txt
fnames=${analysisdir}/bpfnames.txt
bpinfo=${analysisdir}/bpinfo.txt

rm ${tnames}
samtools view $bam NW_003614672.1:181911-183911 | \
  awk '{ print "Insertion1\t" $1 }' >> ${tnames}
samtools view $bam NW_003615350.1:171371-173371 | \
  awk '{ print "Insertion2\t"$1 }' >> ${tnames}
samtools view $bam NW_003617513.1:18412-20412 | \
  awk '{ print "Insertion3\t"$1 }' >> ${tnames}

awk '{ print $2 }' $tnames> $bpnames

allfq=/atium/Data/Nanopore/Analysis/*cho*/*hits.fastq.gz

cat $allfq | seqtk subseq - ${bpnames} | awk '{ print $2 }' | grep . > $fnames

paste -d "\t" $tnames $fnames > $bpinfo

