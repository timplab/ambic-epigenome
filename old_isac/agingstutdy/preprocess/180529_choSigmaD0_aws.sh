#!/bin/bash
env=aws
tardir=/shared/data/tar/fast5
rawroot=/shared/data/raw
outroot=/shared/data/analysis
metadata=$outroot/choSigmaAgingStudy_data.csv
bcall=`readlink -f ../../../oxford/slurm/bcallWrapper.sh`  
align=../../../oxford/slurm/alignWrapper.sh
nwrapper=`readlink -f ../../../nanopolish/slurm/nanopolishWrapper.sh`
ref=/shared/data/Reference/Cricetulus_griseus_chok1gshd.sigmaIgG.fa
nanopolish=/shared/Code/nanopolish/nanopolish

awk 'NR>1' $metadata | while IFS=$',' read -r -a line
do
  n=$(($n+1)) 
#  if [ $n -ge 10 ];then
#    exit
#  fi
  samp=${line[1]}
  echo $samp
  day=${line[3]}
  lab=${line[2]}
  rep=${line[4]}
  tarflag=${line[7]}
  callflag=${line[8]}
  catflag=${line[9]}
  alignflag=${line[10]}
  indexflag=${line[11]}
  cpgflag=${line[12]}
  gpcflag=${line[13]}
  cpgmbedflag=${line[14]}
  gpcmbedflag=${line[15]}
  if [ "$1" == "test" ];then
    echo ${!2}
    continue
  fi
  if [ "$tarflag" == "n" ];then
    echo "untarring $samp"
    com="$bcall untar -e $env --tar $tardir -i $rawroot -o $outroot -b $samp"
    eval $com
    continue
  fi
  if [ "$callflag" == "n" ];then
    echo "bascalling $samp"
    $bcall basecall -e $env -i $rawroot -o $outroot -b $samp
    continue
  fi
  if [ "$catflag" == "n" ];then
    echo "concatenating $samp fastq files"
    $bcall fqcat -e $env -i $rawroot -o $outroot -b $samp --all
    continue
  fi
  if [ "$alignflag" == "n" ];then
    echo "aligning $samp"
    $align -e $env -d $outroot -b $samp -a ngmlr -r $ref
  fi
  if [ "$indexflag" == "n" ];then
    echo "indexing $samp"
    $nwrapper index -e $env -i $rawroot -o $outroot -b $samp -n $nanopolish
    continue
  fi
  if [ "$cpgflag" == "n" ];then
    echo "methyation calling $samp"
    $nwrapper mcall -e $env -o $outroot -b $samp -n $nanopolish -g $ref -m cpg
  fi
#  if [ "$gpcflag" == "n" ];then
#    echo "methyation calling $samp"
#    $nwrapper mcall -e $env -o $outroot -b $samp -n $nanopolish -g $ref -m gpc
#    continue
#  fi
  if [ "$cpgmbedflag" == "n" ];then
    echo "generating cpg methylation bed files for $samp"
    $nwrapper mbed -e $env -o $outroot -b $samp -m cpg
  fi
  if [ "$gpcmbedflag" == "n" ];then
    echo "generating gpc methylation bed files for $samp"
    $nwrapper mbed -e $env -o $outroot -b $samp -m gpc
  fi


done


