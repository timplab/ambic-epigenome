#!/bin/bash

#argument parsing
while :
do
  case "$1" in 
    -i | --input) #path of the raw fast5 folders
      rawpath=`readlink -f $2`
      shift 2
      ;;
    -o | --outdir) # output dir
      outdir=$2
      shift 2
      ;;
    -f | --flowcell) #flowcell
      flow=$2
      shift 2
      ;;
    -k | --kit) #kit
      kit=$2
      shift 2
      ;;
    -b | --barcoding) #barcode?
      arg="--albacore-argument --barcoding"
      shift 1
      ;;
    -s | --src) #source code path
      srcpath=$2
      shift 2
      ;;
    --listpath) ##optional
      listpath=$2
      shift 2
      ;;
    --marcc) ##marcc option
      slurmopt="--partition=shared --time=4:0:0"
      slurmopt="${slurmopt} --mail-type=end --mail-user=ilee29@jhmi.edu"
      slurmopt="${slurmopt} -c 24"
      slurmtype="--marcc"
      shift 1
      ;;
    --aws) ##aws option
      slurmopt="-c 36"
      slurmtype="--aws"
      shift 1
      ;;
    *) break
      ;;
  esac
done

if [ -z $flow ];then
  echo "no flowcell provided, using FLO-MIN106"
  flow=FLO-MIN106
fi
if [ -z $kit ]; then
  echo "no kit provided, using SQK-LSK108"
  kit=SQK-LSK108
fi
if [ ! -d $outdir ]; then
  mkdir $outdir
fi

if [ $listpath ];then
  echo "path list provided in $listpath"
  arraylist=`cat $listpath`
else
  arraylist=`printf "%s\n" $rawpath/*`
fi
n=`echo $arraylist | wc -w`
array="1-$n"

#check
#echo "$arraylist"
#echo $rawpath

echo "following command is used:"
echo "sbatch -a $array -D $outdir $slurmopt \
  $srcpath/oxford_call.scr -i ${rawpath} -o $outdir \
  -f $flow -k $kit --arrayval "$arraylist" $arg $slurmtype"
sbatch -a $array -D $outdir $slurmopt \
  $srcpath/oxford_call.scr -i ${rawpath} -o $outdir \
  -f $flow -k $kit --arrayval "$arraylist" $arg $slurmtype
  
