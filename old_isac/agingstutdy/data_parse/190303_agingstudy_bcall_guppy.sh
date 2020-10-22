#!/bin/bash
rawroot=/data/ambic
broot=/data/basecalled
cd $rawroot
dirs=$(find * -maxdepth 0 -type d)
cfg=/home/grid/Code/ilee/oxford/config/guppy_r941_10kchunk.cfg
dev=1
for dir in $dirs; do
  bdir=$broot/$dir
  log=$broot/$dir.bcall.log
  com="guppy_basecaller --verbose_logs -x cuda:$dev \
    -i $dir -s $bdir --recursive -q 0 -c $cfg &> $log"
  if [[ $dev -eq 1 ]];then
    coms1="$com; $coms1"
    dev=0
  elif [[ $dev -eq 0 ]];then
    coms0="$com; $coms0"
    dev=1
  fi
done
eval $coms1 & PID1=$!
eval $coms0 & PID0=$!
wait $PID1
wait $PID0
