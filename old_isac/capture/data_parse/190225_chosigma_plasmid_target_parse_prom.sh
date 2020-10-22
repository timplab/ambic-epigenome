#!/bin/bash
root=/data/ambic
rdir=$root/reads
samp=190225_choSigma_plasmid_target
bdir=$root/bcall/$samp
[ -e $bdir ]||mkdir -p $bdir
ddir=$root/demux/$samp
[ -e $ddir ]||mkdir -p $ddir

sdir=$rdir/$samp
f5dir=$root/reads_split/$samp
f5demux=$root/reads_demux/$samp
[ -e $f5demux ]||mkdir -p $f5demux

if [ "$1" == "bcall" ];then
  log=$bdir/$samp.bcall.log
  config=/home/prom/Code/ilee/oxford/config/guppy_r941_10kchunk.cfg
  guppy_basecaller -q 0 --num_callers 16 --device cuda:2 \
    -i $sdir -s $bdir -c $config &> $log
fi

if [ "$1" == "demux" ];then
  log=$ddir/$samp.demux.log
  guppy_barcoder -i $bdir -s $ddir -t 90 &> $log
fi

if [ "$1" == "splitf5" ];then
  multi_to_single_fast5 --input_path $sdir --save_path $root/reads_split/$samp
fi

if [ "$1" == "mergefq" ];then
  for dir in $(find $ddir -type d -name "barcode*"); do
    bc=$(basename "$dir")
    fqs=$(find $dir/* -name "*fastq")
    cat $fqs > $ddir/$bc.fastq
  done
fi

if [ "$1" == "index" ];then
  cd $ddir
  for fq in $(find $ddir -maxdepth 1 -type f -name "*fastq"); do
    nanopolish index -d ../../reads_split/190225_choSigma_plasmid_target/ $fq
  done
fi

if [ "$1" == "f5classify" ];then
  scr=/home/prom/Code/ilee/oxford/fast5_split_barcodes.py
  python3 $scr -v -b $ddir/barcoding_summary.txt \
    -r $ddir/barcode01.fastq.index.readdb -n 4000 -o $f5demux
fi
  
fdir=$root/reads_final/$samp
[ -e $fdir ]||mkdir $fdir
if [ "$1" == "mergef5" ];then
  for dir in $(find $f5demux/* -maxdepth 0 -type d); do
    bc=$(basename $dir)
    single_to_multi_fast5 -t 36 --recursive -i $dir -s $fdir/$bc -n 10000
  done
fi
