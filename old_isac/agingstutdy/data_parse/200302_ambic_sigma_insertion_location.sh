#!/bin/bash
assembly=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.fa
mmi=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.mapont.mmi

ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
refmmi=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.mmi

root=/kyber/Data/Nanopore/projects/ambic/sigma
fqdir=$root/pooled_reps/fastq
bamdir=$root/bam_plasmid

# first subset reads that map to the plasmid at all
if [ "$1" == "subset" ];then
  for fq in $(find $fqdir -name "*fastq.gz" ); do
    base=$(basename "${fq%%.*}")
    echo $base
    bam=$bamdir/$base.plasmid.bam
    fqsub=$bamdir/$base.plasmid.fastq
    log=$bamdir/$base.plasmid.log
    minimap2 -t 68 -a $mmi $fq 2> $log |\
      samtools view -@ 12 -F 4 -hb - > $bam
    samtools fastq $bam > $fqsub
  done
fi

# then let's map these to the plasmid again, outputting paf
if [ "$1" == "map" ];then
  for fq in $(find $bamdir -name "*fastq" ); do
    base=$(basename "${fq%%.*}")
    bases="$bases $base"
  done
  base="{}"
  fq=$bamdir/$base.plasmid.fastq
  plasmidpaf=$bamdir/$base.plasmid.paf
  com="minimap2 -c $mmi $fq > $plasmidpaf"
  parallel $com ::: $bases
  genomepaf=$bamdir/$base.genome.paf
  com="minimap2 -c $refmmi $fq > $genomepaf"
  parallel $com ::: $bases
fi

# parse pafs
if [ "$1" == "parse" ];then
  for fq in $(find $bamdir -name "*fastq" ); do
    base=$(basename "${fq%%.*}")
    echo $base
    plasmidpaf=$bamdir/$base.plasmid.paf
    genomepaf=$bamdir/$base.genome.paf
    out=$bamdir/$base.pafparse.out
    python ./parse_insertion_paf.py $plasmidpaf $genomepaf > $out 

  done
fi



