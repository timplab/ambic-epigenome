#!/bin/bash
root=/shared
bamroot=${root}/bam
wdir=${root}/insertion
aligner=ngmlr
########aligner=minimap2
bamdir=${bamroot}/$aligner
svroot=${root}/sniffles
svdir=${svroot}/plasmid
[ -e $svdir ]||mkdir -p $svdir
ref=/shared/ref/criGri1_plasmid2.fa
slurmscr=/shared/Code/ilee/sv/slurm/sniffles.scr
srcdir=/shared/Code
pname=4_0cdhfr_vrc01wtg1m3_dgv 
insertdir=${root}/insertion
inbamchunks=$insertdir/bam_chunks
inbamdir=$insertdir/bam
rnamedir=$insertdir/readnames
fqdir=$insertdir/fastq
outbamdir=$insertdir/insertbam

## extracting inserrtion reads
#for bam in `ls $bamdir/*bam`;do
#  base=$(basename "$bam")
#  base=${base%.sorted.bam}
#  samp=${base%%.*}
#  ind=${base#*.}
#  out=$inbamchunks/$base.insertion.bam
#  echo $base
#  samtools view -F 2048 -hb $bam $pname > $out
#done

## merge same sample
#for dir in `printf "%s\n" /shared/raw/*`;do
#  samp=$(basename "$dir")
#  bams=`find $inbamchunks -name "$samp*bam"`
#  echo $samp
#  samtools merge -f $inbamdir/$samp.insertion.bam $bams
#done

### redo mapping on merged samples
#align=/shared/Code/ilee/oxford/slurm/oxford_align.scr
#aligner=ngmlr
#ref=/shared/ref/picr_IgG2.fa
#for bam in `ls $inbamdir/*bam`
#do
#  base=$(basename "$bam")
#  base=${base%.insertion.bam}
#  echo $base
#  fq=$fqdir/$base.insertion.fq
#  prefix=$outbamdir/$base.insertion
#  log=$prefix.align.log
#  echo $log,$fq,$prefix
#  bedtools bamtofastq -i $bam -fq $fq
#  sbatch -o $log -e $log \
#    $align -i $fq -b $prefix -a $aligner -r $ref --aws
#done

## extract read names and fastq seq of insertion reads and redo mapping
#align=/shared/Code/ilee/oxford/slurm/oxford_align.scr
#aligner=ngmlr
#ref=/shared/ref/picr_IgG2.fa
#for bam in $inbamchunks/choNIHhost1.113*bam #`ls $inbamchunks/*bam`
#do
#  base=$(basename "$bam")
#  base=${base%.insertion.bam}
#  echo $base
#  fq=$fqdir/$base.insertion.fq
##  samtools view $bam | cut -f1 | uniq > $rnamedir/$base.readnames.txt
##  bedtools bamtofastq -i $bam -fq $fq
#  prefix=$outbamdir/$base.insertion
#  log=$prefix.align.log
#  sbatch -o $log -e $log \
#    $align -i $fq -b $prefix -a $aligner -r $ref --aws
#done

## indexing
#dbdir=/shared/fastq
#logdir=/shared/log/insertindex
#polpath=/shared/Code/nanopolish/nanopolish
#rawroot=/shared/raw
#for fq in `find $fqdir -name "*fq"`
#do
#  base=$(basename "$fq")
#  samp=${base%.insertion.fq}
#  bam=`readlink -f $inbamdir/$samp*`
#  rep=${samp%%.*}
#  ind=${samp##*.}
#  rawdir=$rawroot/$rep/$ind
#  logpath=$logdir/$samp.index.log
#  #echo $fq #$rawdir,$logpath,$fq
#  ## first step is used only to get the fai stuff
##  sbatch -o $logpath -e $logpath $srcdir/ilee/dnamods/slurm/fast5index.scr -s $polpath -i $fq -d /shared/tmp
#  ## subset thre readdb 
#  samtools view $bam | cut -f1 > $rnamedir/$samp.readnames.txt
#  rname=`readlink -f $rnamedir/$samp.readnames.txt`
#  db=`readlink -f /shared/fastq/$samp.readdb`
#  /shared/Code/ilee/nanopolish/subsetindex.py -r $db -n $rname > $wdir/readdb/$base.readdb
#done
#mv $wdir/readdb/* $wdir/fastq
#
## methylation calling
root=/shared
wdir=/shared/insertion
logdir=/shared/log/insertmeth
[ -e $logdir ]||mkdir $logdir
methdir=$wdir/meth
[ -e $methdir ]||mkdir $methdir
srcpath=$root/Code/ilee/dnamods
poldir=$root/Code/nanopolish
polpath=$root/Code/nanopolish/nanopolish
ref=$root/ref/picr_IgG2.fa
for fq in `find $fqdir -name "*fq"`;do
  base=$(basename "$fq")
  samp=${base%.insertion.fq}
  bam=`readlink -f $outbamdir/$samp*sorted.bam`
  logpath=$logdir/$samp.methcall.log
  methout=$methdir/$samp.insertion.meth.tsv
  sbatch -o $logpath -e $logpath $srcpath/slurm/callmeth.scr \
    -s $polpath -r $fq -g $ref -b $bam -o $methout
done


