#!/bin/bash
wdir=/atium/Data/Nanopore/Analysis/170608_cho
rawdir=/atium/Data/Nanopore/oxford
base=choIgGNIH
fadir=${wdir}/fasta
blastdir=${wdir}/blast
csvpath=${wdir}/*${base}_all.csv.gz
outpath=${wdir}/hitfast5
tmpdir=~/tmp/cho
mkdir $tmpdir

for blastpath in `readlink -f $blastdir/*plasmid_v2*`
do
  blast=$(basename $blastpath)
  samp=${blast%%_plasmid*}
  echo $samp
  fapath=`readlink -f $fadir/$samp.fasta`

  echo "get template name"
  awk '{ print $1 }' $blastdir/$blast | uniq > $tmpdir/${samp}_tnames.txt

  echo "get file name"
  seqtk subseq $fapath $tmpdir/${samp}_tnames.txt | \
    awk '{ print $2}' | \
    grep . > ${tmpdir}/${samp}_fnames.txt
  
  echo "get fast5 files"
  tarpath=`find -L ${rawdir}/${samp}/ -type f \( -iname "${samp}*tar.gz" ! -iname "*log*" \)`
  for t in $tarpath; do
    pathtar=${tmpdir}/${samp}_tarpath.txt
    fnames=${tmpdir}/${samp}_fnames.txt
    tar -tzf $t | grep -f $fnames > $pathtar
    echo "extracting $t"
    tar -xzf $t -C $outpath --files-from $pathtar
  done
done

#rm -r $tmpdir
