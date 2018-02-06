#!/bin/bash
while :
do
  case "$1" in
    -s | --src) #path of source code
      srcpath=$2
      shift 2
      ;;
    -w | --wdir) #path of working dir
      wdir=$2
      shift 2
      ;;
    -a | --array) #subset 
      array=$2
      shift 2
      ;;
    --listpath)#optional
      listpath=$2
      shift 2
      ;;
    -b | --barcoding) #barcode
      extarg="$extarg $1"
      shift 1
      ;;
    --marcc) #marcc specific options
      extarg="$extarg $1"
      slurmopt="--partition=shared --time=30:0"
      slurmopt="$slurmopt --mail-type=end --mail-user=ilee29@jhmi.edu"
      slurmtype="--marcc"
      shift 1
      ;;
    *) break
      ;;
  esac
done

callpath=${wdir}/called
ext=${wdir}/ext

if [ ! -d $ext ];then
  mkdir $ext
fi

if [ $listpath ]; then
  echo "path list provided in $listpath"
  arraylist=`cat $listpath`
else
  pwd=$PWD
  cd $callpath
  arraylist=`ls -d */`
  cd $PWD
fi
n=`echo $arraylist | wc -w`
array="1-$n"

##check
#echo $extarg
#echo $array
echo "The command to be performed:"
echo "sbatch -a $array -D $ext $slurmopt \
  $srcpath/extract_call.scr -s $srcpath -w $wdir --arrayval "$arraylist" $extarg"
sbatch -a $array -D $ext $sbatchopt \
  $srcpath/extract_call.scr -s $srcpath -w $wdir --arrayval "$arraylist" $extarg
