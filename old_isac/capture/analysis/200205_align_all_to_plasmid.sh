#!/bin/bash
assembly=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_singler_2line.fa
mmi=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_singler_2line.mapont.mmi
ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
refmmi=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.mmi
root=/kyber/Data/Nanopore/projects/ambic/capture/190225_choSigma_plasmid_target_pipeline_test/data/nanopore
fqdir=$root/reads
outdir=$root/align_plasmid

if [ ! -e $mmi ];then
  minimap2 -x map-ont -d $mmi $assembly

fi


for fq in $(find $fqdir -type f -name "*fastq.gz"); do
  base=$(basename "${fq%.fastq.gz}")
  out=$outdir/$base.plasmid.bam
  log=$outdir/$base.plasmid.log
  # align
#  minimap2 --MD -L -t 32 -a $mmi $fq 2> $log |\
#    samtools view -b - | samtools sort -o $out
#  samtools index $out
#  samtools idxstats $out > $outdir/$base.plasmid.idxstats.txt

  # extract fastq of plasmid-aligned data
  outfq=$outdir/$base.plasmid.fastq
#  samtools view $out polished_singler -b | samtools fastq - > $outfq

  # align to the genome
  log=$outdir/$base.plasmid.ref.log
  refout=$outdir/$base.plasmid.ref.bam
#  minimap2 --MD -L -a $refmmi $outfq 2> $log |\
#    samtools view -b - | samtools sort -o $refout
#  samtools index $refout
  
  # then get out the coordinates of alignments
  bed=$outdir/$base.plasmid.ref.bed
#  bedtools bamtobed -i $refout > $bed
done

# process the outputs
Rscript 200205_targeted_insertion_candidates.Rmd

