#!/bin/bash
assembly=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.fa
mmi=/kyber/Data/Nanopore/projects/ambic/sigma/plasmid_assembly/polished_tripler_twoline.mapont.mmi

ref=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.fa
refmmi=/mithril/Data/NGS/Reference/cho/picr_ensembl/Cricetulus_griseus_picr.CriGri-PICR.dna.toplevel.mmi

root=/kyber/Data/Nanopore/projects/ambic/capture/190225_choSigma_plasmid_target_pipeline_test/data/nanopore
fqdir=$root/reads
outdir=$root/align_plasmid

if [ ! -e $mmi ];then
  minimap2 -x map-ont -d $mmi $assembly
  minimap2 -x map-ont -d $refmmi $ref
fi


if [ "$1" == "align" ];then
  for fq in $(find $fqdir -type f -name "*fastq.gz"); do
    base=$(basename "${fq%.fastq.gz}")
    echo $base
    plasmidout=$outdir/$base.plasmid.paf
    refout=$outdir/$base.reference.paf
    log=$outdir/$base.plasmid.log
    # get paf to the plasmid
    minimap2 -x ava-ont $assembly $fq 2> $log > $plasmidout
    # also get paf to reference
    log=$outdir/$base.reference.log
    minimap2 -t 24 -c $refmmi $fq 2> $log > $refout
  done
fi



if [ "$1" == "flank" ];then
  flank=$outdir/200227_ambic_targeted_genomic_alignments_genomic_flanks.txt
  [ ! -e $flank ]||rm $flank
  for fq in $(find $fqdir -type f -name "*fastq.gz"); do
    base=$(basename "${fq%.fastq.gz}")
    plasmidout=$outdir/$base.plasmid.paf
    refout=$outdir/$base.reference.paf
    # let's find reads that have flanking sequence that does not align to assembly at all
    python ./find_flanking_from_paf.py $plasmidout $refout |\
      awk -v samp=$base 'OFS = "\t"{ print $0,samp }' >> $flank
  done
  # use these flanking alignments to get candidate regions
  # 1kb region outside of the candidate insertion site for idt
  grna=$outdir/200227_ambic_gRNA_coordinates.txt
  python ./get_candidate_coord_from_flank_regions.py $flank  > $grna
  # let's get sequences in these regs and input them into idt
fi

if [ "$1" == "seq" ];then
  grna=$outdir/200227_ambic_gRNA_coordinates.txt
  seq=$outdir/200227_ambic_gRNA_sequences.fa
  [ ! -e $seq ]||rm $seq
  while IFS=$'\t' read -r -a line || [[ -n "$line" ]]
  do
    coord=${line[5]}
    samtools faidx $ref $coord >> $seq
  done <<< "$(awk 'NR>1' $grna)"
  sed -i -e "s|:|_|g" $seq
fi


