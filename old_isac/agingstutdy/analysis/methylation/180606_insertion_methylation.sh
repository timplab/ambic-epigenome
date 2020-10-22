#!/bin/bash
# find methylation at and/or around insertion sites
day=0
root=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/D$day
methdir=$root/methylation
poolmethdir=$root/pooled/methylation/methbyread

region="scaffold_16:15304132-15305120"
startcoord=15304132
stopcoord=15305120
reglab="scaf16"
reg=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/insertion/scaffold_16_insertion.bed

#region="ambic_sigma_IgG_HC"
#startcoord=1
#stopcoord=1500
#reglab="sigmaIgGHC"
#reg=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/insertion/sigmaIgGHC.bed

parser=/home/isac/Code/ilee/nanopolish/script/parseMethylbed.py
plotter=/home/isac/Code/ilee/plot/methylation_average_plots.R
plotdir=/home/isac/Dropbox/Data/ambic/aging_study/plots

insertdir=$root/insertion

if [ "$1" == "qc" ];then
  freqs=`find $methdir -name "*$reglab*freq.txt"`
  for freq in $freqs;do
    echo $freq
    awk '{ meth+=$4;unmeth+=$5 }END{ print meth,unmeth,meth/(meth+unmeth) }' $freq
  done
  exit
fi
if [ "$1" == "plot" ];then
  freqs=`find $methdir -name "*${reglab}*insert.methfreq.txt"`
  echo $freqs
  out=$plotdir/$reglab.methfreq.pdf
  Rscript $plotter methByRegion -o $out -r $reg --coverage $freqs 
  exit
fi

for samp in Host StableGlut UnstableGlut;do
  lab=${samp}D$day
  if [ "$1" == "insertion" ];then
  # first let's get methylation of the insertion-containing reads
    bed=`find $insertdir -name "$lab*bed.gz"`
    out=$methdir/$lab.$reglab.insert.methfreq.txt
    tabix $bed $region | cut -f4 | uniq | wc -l #python $parser frequency > $out
  elif [ "$1" == "not" ];then
  # this region in reads not in insertion
    insertnames=$insertdir/$lab.insertion.rnames.txt
    wc -l $insertnames
    beds=`find $poolmethdir -name "$lab*bed.gz"`
    for bed in $beds;do
      base=$(basename "$bed")
      pre=${base%%.*}
      echo $pre
      regbed=$methdir/$pre.$reglab.methbyread.bed
      cut -f4 $regbed | sort | uniq | wc -l
      tabix $bed $region > $regbed
    done
    out=$methdir/$lab.$reglab.notinsert.methfreq.txt
    cat $methdir/$lab*methbyread.bed |\
      grep -F -v -f $insertnames |\
      sort -k1,1 -k2,2n |\
      python $parser frequency > $out
  fi
done

