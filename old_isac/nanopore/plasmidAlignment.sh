#!/bin/bash
root=/atium/Data/Nanopore/Analysis
rawdir=/atium/Data/Nanopore/oxford
refpath=/mithril/Data/NGS/Reference/cho/criGri1_IgG.fa
allpath=/atium/homes/isac/Data/ambic/nanopore
allbam=${allpath}/170608_choIgGNIH_hits.bam
sortbam=${allpath}/170608_choIgGNIH_hits.sorted.bam

for runpath in `ls -d ${root}/*cho*`; do
  base=$(basename "$runpath")
  echo $base
  hitfastq=`ls ${runpath}/${base}_hits.fastq.gz`
  hitsam=${runpath}/${base}_hits.sam
  hitbam=${runpath}/${base}_hits.bam

  bwa mem -x ont2d $refpath $hitfastq > ${runpath}/${base}_hits.sam
  samtools view -b $hitsam -o $hitbam
  rm $hitsam

  echo "done with sample $base"

done

echo "concatenating the bam"
samtools cat -o $allbam ${root}/*cho*/*bam
cp ${root}/*cho*/*bam $allpath

samtools sort -o $sortbam $allbam
rm $allbam
samtools index $sortbam
