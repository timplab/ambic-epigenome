#!/bin/bash
f5root=/data/ambic
bcalldir=/data/basecalled
for dir in $(find $f5root/* -maxdepth 0 -type d); do
  base=$(basename "$dir")
  bases="$bases $base"
done
if [ "$1" == "bcall" ];then
  devs=$(seq 0 3)
  for dev in $devs; do
    script=$bcalldir/bcall_cuda$dev.sh
    echo "#!/bin/bash" > $script
    chmod a+x $script
  done
  i=0
  for base in $bases; do
    outdir=$bcalldir/$base
    if [ -e $outdir ];then
      continue
    fi
    echo $base
    f5s=$(find $f5root/$base -type f -name "*fast5")
    if [[ $f5s == *"PAD"* ]]; then
      cfg=/home/prom/Code/ilee/oxford/config/guppy_r941_10kchunk_prom.cfg
    else
      cfg=/home/prom/Code/ilee/oxford/config/guppy_r941_10kchunk.cfg
    fi
#    args="$args $base $cfg $i"
    indir=$f5root/$base
    log=$bcalldir/$base.bcall.log
    com="guppy_basecaller --verbose_logs -x cuda:$i \
      -i $indir --num_callers 16  -s $outdir --recursive -q 0 -c $cfg &> $log"
    script=$bcalldir/bcall_cuda$i.sh
    echo $com >> $script
    if [ $i -eq 3 ];then
      i=0
    else
      i=$(($i+1))
    fi
  done
  com="$bcalldir/bcall_cuda{}.sh"
  parallel $com ::: $devs
fi

if [ "$1" == "cat" ];then
  outdir=/data/sorted
  for base in $bases; do
    echo $base
    dir=$bcalldir/$base
    fqs=$(find $dir -name "*fastq")
    sum=$(find $dir -name "sequencing_summary.txt")
    outfq=$outdir/$base.fastq.gz
    outsum=$outdir/$base.summary.txt
    cp $sum $outsum
    cat $fqs | gzip > $outfq
#    awk 'NR==1||$1!="filename"{ print}' $sums > $outsum
  done
fi

