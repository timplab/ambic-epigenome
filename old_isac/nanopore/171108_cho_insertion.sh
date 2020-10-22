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

for bam in `ls $bamdir/*bam`;do
  echo $bam
done

##getting fastq sequences of the reads that overlap the insertion poitns
#for bam in `ls $bamdir/*pooled.bam`
#do
#  samp=$(basename "$bam")
#  samp=${samp%%.*}
#  echo $samp
#  samtools view -F 2048 -bh $bam \
#    4_0cdhfr_vrc01wtg1m3_dgv \
#    KE382060:5042401-5052401 \
#    KE382060:6621603-6631603 \
#    KE378516:15159-25159 \
#    KE379390:121159-131159 \
#    > $wdir/$samp.insertion.bam
#  ins=$wdir/$samp.insertion.bam
#  bedtools bamtofastq -i $ins -fq $wdir/$samp.insertion.fq
#done

## perform fast5 indexing
#logdir=/shared/log
#srcpath=/shared/Code/ilee/dnamods
#rawroot=/shared/raw
#polpath=/shared/Code/nanopolish/nanopolish
#
## this fails because of subseeting - I use this to create the indexes and for readdb subsetting, go below.
#for samp in host IgG
#do
#  fq=`ls $wdir/*$samp*insertion.fq`
#  base=$(basename "$fq")
#  base=${base%%.*}
#  rawdir=$rawroot/$samp
#  logpath=$logdir/${base}.fast5index.log
#  echo $fq
#  echo "sbatch -o $logpath -e $logpath $srcpath/slurm/fast5index.scr -s $polpath -i $fq -d $rawdir"
#  sbatch -o $logpath -e $logpath $srcpath/slurm/fast5index.scr -s $polpath -i $fq -d $rawdir
#done
#

#db=/shared/fastq
### get readnames and extract readdb read loci
#for bam in `ls $wdir/*bam`
#do
#  base=$(basename "$bam")
#  base=${base%%.*}
#  rn=$wdir/$base.insertion.readnames
#  echo $base
##  samtools view -F 2048 $bam | awk '{ print $1 }' > $rn
#  echo "readdb $base" 
##  cat $db/$base*readdb | grep -f $rn > $wdir/$base.fq.index.readdb
##  /shared/Code/ilee/nanopolish/subsetindex.py -r $db/$base.pooled.readdb \
##    -n $rn > $wdir/$base.insertion.fq.index.readdb2
#done
#
###let's try methylation calling
#logdir=$root/log
#methdir=$root/meth
#srcpath=$root/Code/ilee/dnamods
#poldir=$root/Code/nanopolish
#polpath=$root/Code/nanopolish/nanopolish
#ref=$root/ref/criGri1_plasmid2.fa
#for samp in IgG
#do
#  fq=`ls $wdir/*$samp*.fq`
#  bam=`ls $wdir/*$samp*bam`
#  base=$(basename "$fq")
#  base=${base%%.*}
#  echo $base
#  logpath=$logdir/$base.insertion.methcall.log
#  methout=$methdir/$base.insertion.meth.tsv
#  echo "sbatch -o $logpath -e $logpath $srcpath/slurm/callmeth.scr \
#    -s $polpath -r $fq -g $ref -b $bam -o $methout"
##  sbatch -o $logpath -e $logpath $srcpath/slurm/callmeth.scr \
##    -s $polpath -r $fq -g $ref -b $bam -o $methout
#done
#
## methylation frequency
for samp in host IgG
do
  fq=`ls $wdir/*$samp*.fq`
  bam=`ls $wdir/*$samp*bam`
  base=$(basename "$fq")
  base=${base%%.*}
  echo $base
  methout=$methdir/$base.insertion.meth.tsv
  python $poldir/scripts/calculate_methylation_frequency.py \
    -c 2.5 -i $methout > $methdir/$base.insertion.methfreq.tsv
done


