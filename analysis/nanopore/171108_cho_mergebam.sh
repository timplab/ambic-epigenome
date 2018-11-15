#SBATCH
#SBATCH --job-name=mergebam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-puer-task=16

root=/shared
bamroot=$root/bam
aligner=ngmlr
bamdir=$bamroot/$aligner
bamout=$root/mergedbam
src=/shared/Code
logdir=/shared/log/mergebam
[ -e $logdir ]||mkdir $logdir
mergescr=$src/ilee/slurm/mergebam.scr

#for dir in `find /shared/raw/* -maxdepth 0 -type d`
#do
#  samp=$(basename "$dir")
#  echo $samp
#  p=$bamdir/*$samp*sorted.bam
#  out=$bamout/choNIH$samp.pooled.bam
#  logpath=${logdir}/$samp.merge.log
#  bams=`ls $p`
#  echo "sbatch --output=$logpath --error=$logpath $mergescr -s $src -o $out -p "$p""
#  sbatch --output=$logpath --error=$logpath $mergescr -s $src -o $out -p "$p"
#done
## second round of merging
for samp in host IgG
do
  echo $samp
  p=$bamout/*$samp*sorted.bam
  out=$bamout/choNIH$samp.pooled.bam
  logpath=${logdir}/$samp.merge.log
  bams=`ls $p`
  echo "sbatch --output=$logpath --error=$logpath $mergescr -s $src -o $out -p "$p""
  sbatch --output=$logpath --error=$logpath $mergescr -s $src -o $out -p "$p"
done

