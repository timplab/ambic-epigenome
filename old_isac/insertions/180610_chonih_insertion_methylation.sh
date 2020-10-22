#!/bin/bash
# find methylation at and/or around insertion sites
root=/dilithium/Data/Nanopore/Analysis/171025_cho
methdir=$root/meth/readlevel

region="picr_49:6218858-6220858"
startcoord=6218858
stopcoord=6220858
reglab="picr_49"
reg=$root/insertion/picr/picr49_insertion.bed

#region="ambic_sigma_IgG_HC"
#startcoord=1
#stopcoord=1500
#reglab="sigmaIgGHC"
#reg=/dilithium/Data/Nanopore/projects/choSigmaAgingStudy/analysis/insertion/sigmaIgGHC.bed

parser=/home/isac/Code/ilee/nanopolish/script/parseMethylbed.py
plotter=/home/isac/Code/ilee/plot/methylation_average_plots.R
plotdir=/home/isac/Dropbox/Data/ambic/plots

insertdir=$root/insertion/picr

if [ "$1" == "qc" ];then
  freqs=`find $methdir -name "*$reglab*freq.txt"`
  for freq in $freqs;do
    echo $freq
    awk '{ meth+=$4;unmeth+=$5 }END{ print meth,unmeth,meth/(meth+unmeth) }' $freq
  done
  exit
fi
if [ "$1" == "plot" ];then
  freqs=`find $insertdir -name "*${reglab}*insert.methfreq.txt"`
  echo $freqs
  out=$plotdir/$reglab.methfreq.pdf
  Rscript $plotter methByRegion -o $out -r $reg --coverage $freqs 
  exit
fi

for samp in host IgG;do
  lab=choNIH$samp
  if [ "$1" == "insertion" ];then
  # first let's get methylation of the insertion-containing reads
    bed=`find $insertdir -name "$lab*bed.gz"`
    echo $bed
    out=$insertdir/$lab.$reglab.insert.methfreq.txt
    tabix $bed $region | python $parser frequency > $out
  elif [ "$1" == "not" ];then
  # this region in reads not in insertion
    insertnames=$insertdir/picr49.insertion.rnames.txt
    wc -l $insertnames
    beds=`find $methdir -name "$lab*bed.gz"`
    echo $beds
    for bed in $beds;do
      base=$(basename "$bed")
      pre=${base%%.*}
      echo $pre
      regbed=$insertdir/$pre.$reglab.methbyread.bed
      tabix $bed $region > $regbed
      cut -f4 $regbed | sort | uniq | wc -l
    done
    out=$insertdir/$lab.$reglab.notinsert.methfreq.txt
    cat $insertdir/$lab*methbyread.bed |\
      grep -F -v -f $insertnames |\
      sort -k1,1 -k2,2n |\
      python $parser frequency > $out
  fi
done

