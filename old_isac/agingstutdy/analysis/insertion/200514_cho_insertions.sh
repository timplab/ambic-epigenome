#!/bin/bash
day=0
root=/kyber/Data/Nanopore/projects/ambic/sigma
bamdir=$root/bam/pooled_rep
outdir=$root/insertion
insertiondir=$root/insertion
[ -e $insertiondir ]||mkdir $insertiondir

assembly=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.fa
mmi=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.mapont.mmi

ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
refmmi=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.mmi

for samp in Host StableGln StableNogln UnstableGln UnstableNogln;do
  for day in 0 30 60 90; do
    for rep in 1 2 3; do
        label=CHOZN${samp}Day${day}_$rep
        labs="$labs $label"
      done
    done
done
label="{}"
bam=$bamdir/$label.pooled.bam

# first let's get names of reads that aligned to the plasmid as well as the genome
insertreads=$insertiondir/$label.insertion.bam
insertfq=$insertiondir/$label.insertion.fastq
if [ "$1" == "getreads" ];then
  com="samtools view -h $bam \
  ambic_sigma_IgG_LC ambic_sigma_IgG_HC |\
  samtools view  -hb \
  > $insertreads"
  parallel $com ::: $labs
fi

# number of reads for each sample that aligned to th plasmid?
if [ "$1" == "readqc" ];then
  for lab in $labs; do
    r=$insertiondir/$lab.insertion.bam
    samtools view $r | wc -l | awk -v name=$lab 'OFS="\t"{ print name,$0}' | tr "_" "\t" 
  done
fi

if [ "$1" == "getfq" ];then
  com="samtools fastq $insertreads > $insertfq"
  parallel $com ::: $labs
fi

# get paf
if [ "$1" == "paf" ];then
  fq=$insertfq
  plasmidpaf=$insertiondir/$label.plasmid.paf
  com="minimap2 -c $mmi $fq > $plasmidpaf"
  parallel $com ::: $labs
  genomepaf=$insertiondir/$label.genome.paf
  com="minimap2 -t 8 -c $refmmi $fq > $genomepaf"
  parallel -j 6 $com ::: $labs
fi

if [ "$1" == "flank" ];then
  scr="/home/isac/Code/ilee/project/ambic/capture/analysis/find_flanking_from_paf.py"
  plasmidpaf=$insertiondir/$label.plasmid.paf
  genomepaf=$insertiondir/$label.genome.paf
  com="$scr $plasmidpaf $genomepaf"
  candidates=$insertiondir/insertion_point_candidates.txt
  [ ! -e $candidates ]||rm $candidates

  for lab in $labs; do
    plasmidpaf=$insertiondir/$lab.plasmid.paf
    genomepaf=$insertiondir/$lab.genome.paf
    com="$scr $plasmidpaf $genomepaf"
    $com | awk -v name=$lab 'OFS="\t"{ print name,$0}' | tr "_" "\t" >> $candidates
  done
fi

# let's summarize the insertion points
if [ "$1" == "summarize" ];then
  Rscript ./200611_summarize_insertion.R
fi



