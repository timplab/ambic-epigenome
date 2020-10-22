#!/bin/bash
srcdir=$(dirname "$0")
coderoot=$(readlink -f "$srcdir/../../../../")
export LC_ALL=C
binroot="/mnt/ambic1"
export PATH="$binroot/bin:$binroot/Code/miniconda3/bin:$PATH"
fqdir=$binroot/fastq
root=/shared/Data
calldir=$root/bcall
droot=$binroot/demux
bamdir=$root/bam
mroot=$root/mcall
mbeddir=$root/mbed

samps="$binroot/ambic_samples.csv"
# for bcall, use bash scripts
if [ "$1" == "bcall" ];then
  echo "bcall"
  wrapper="$coderoot/oxford/slurm/call_wrapper.sh"
  awk 'NR>1' $samps | while IFS=$',' read -r -a line || [[ -n "$line" ]]; do
    run=${line[0]}
    dir=$(find /mnt -maxdepth 3 -type d -name "$run")
    fcell=${line[2]}
    if [ "$fcell" == "FLO-PRO002" ];then
      config="$coderoot/oxford/config/guppy_10kchunk_prom.cfg"
    elif [ "$fcell" == "FLO-MIN106" ];then
      config="$coderoot/oxford/config/guppy_10kchunk.cfg"
    fi
    bcallargs="-c $config"
    outdir=$broot/$run
    echo $outdir
#    [ -e $outdir ]||mkdir -p $outdir
    $wrapper -i $dir -o $outdir -a "$bcallargs" -e "aws" --guppy --dryrun
  done
fi

if [ "$1" == "cat" ];then
  echo "cat fastqs"
fi

# demux reads
if [ "$1" == "demux" ];then
  echo "demux"
  awk 'NR>1' $samps | while IFS=$',' read -r -a line || [[ -n "$line" ]]; do
    run=${line[0]}
    barcode=${line[4]}
    if [ "$barcode" == "N" ];then
      continue
    fi
    fq=$(find $fqdir/$run -name "$run.fastq.gz")
    if [ -z $fq ] ;then
      echo $run
      outdir=$droot/$run
      log=$outdir/$run.demux.log
      [ -e $outdir ]||mkdir $outdir
      guppy_barcoder -t 72 -i $fqdir/$run -s $outdir > $log
    fi
  done
fi
