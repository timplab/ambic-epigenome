#!/bin/bash 
root=/shared
codedir=$root/Code
repo=$codedir/ilee
srcpath=$repo/dnamods
slurmpath=$repo/oxford/slurm
polpath=$codedir/nanopolish/nanopolish
npdir=~/Code/nanopolish
logroot=$root/log
[ -e $logdir ]||mkdir $logdir
rawtardir=$root/rawtar
[ -e $rawtardir ]||mkdir $rawtardir
fqdir=$root/fqchunks
refdir=$root/ref
rawdir=$root/raw
bamdir=$root/bam

## first fetch the data first of host
#aws s3 sync s3://timpnanopore/oxford/171025_cho/ $rawtardir --exclude "*" --include "*host*" 

## fetch rest of the data
#aws s3 sync s3://timpnanopore/Analysis/171025_cho/ $root
#aws s3 sync s3://timp2ndgen/Reference/cho/criGri1/ $refdir

### untar raw data
#logdir=$logroot/untar
#[ -e $logdir ]||mkdir -p $logdir
#for tag in IgG host
#do
#  for i in 1 2 3 
#  do
#    repdir=$rawdir/$tag/$i
#    tars=`ls $rawtardir/*$tag*$i*.raw.tgz`
#    for tar in $tars;do
#      base=$(basename "$tar")
#      base=${base%%.*}
#      basedir=$repdir/$base
#      [ -e $basedir ]||mkdir -p $basedir
#      echo $base
#      logpath=${logdir}/${base}.untar.log
#      echo "sbatch -o $logpath -e $logpath -D $basedir $slrumpath/untar.scr $tar"
#      sbatch -o $logpath -e $logpath -D $basedir $slurmpath/untar.scr $tar
#    done
#  done
#done

## perform fast5 indexing
#logdir=$logroot/fast5index
#[ -e $logdir ]||mkdir -p $logdir
#for fq in `ls $fqdir/*.fq.gz`
#do
#  dir=$(dirname "$fq")
#  base=$(basename "$fq")
#  base=${base%.fq.gz}
#  samp=${base%%.*}
#  ind=${base##*.}
#  repdir=$rawdir/$samp/$ind
#  logpath=$logdir/$base.fast5index.log
#  if [ -e  $fq.index.readdb ];then
#    x=x
#    echo "done $base"
#  else
#    echo "sbatch -o $logpath -e $logpath $srcpath/slurm/fast5index.scr -s $polpath -i $fq -d $repdir"
#    sbatch -o $logpath -e $logpath $srcpath/slurm/fast5index.scr -s $polpath -i $fq -d $repdir
#  fi
#done

# perform methylation calling
methdir=$root/meth
[ -e $methdir ]||mkdir $methdir
ref=$refdir/picr_IgG2.fa
logdir=$logroot/methcall
[ -e $logdir ]||mkdir -p $logdir
for db in `ls $fqdir/*.readdb`
do
  base=$(basename "$db")
  base=${base%.fq.gz.index.readdb}
  logpath=$logdir/$base.methcall.log
  fq=`ls $fqdir/$base.fq.gz`
  bam=`ls $bamdir/ngmlr/$base.sorted.bam`
  methout=$methdir/$base.meth.tsv
  echo "sbatch -o $logpath -e $logpath $srcpath/slurm/callmeth.scr \
    -s $polpath -r $fq -g $ref -b $bam -o $methout"
  sbatch -o $logpath -e $logpath $srcpath/slurm/callmeth.scr \
    -s $polpath -r $fq -g $ref -b $bam -o $methout
done

