#!/bin/bash 
## usage : ./untar_subset.sh -i /path/to/tar -n number -o /dir/to/out/ [-s extension]
while :
do
  case "$1" in
    -i | --infile) #path to input tar
      infile=`readlink -f $2`
      shift 2
      ;;
    -n | --numreads) # number of reads to extract
      n=$2
      shift 2
      ;;
    -o | --outdir) # directory of output
      outdir=$2
      shift 2
      ;;
    -s | --suffix ) ## suffix for file
      suffix=$2
      shift 2
      ;;
    -t | --threads ) ##number of instances of tar
      t=$2
      shift 2
      ;;
    -h | --help ) ## display help
      printf "usage: untar_subset.h [-h] --infile /PATH/TO/TGZ --numreads NUM -outdir /DIR/TO/OUT [--suffix EXT ]\
        \n\nExtract a subset of a tar.gz file.\
        \n\narguments:\
        \n  -h, --help\t\t\tshow this help message and exit\
        \n  -i, --infile /PATH/TO/TGZ\tpath to input tar.gz file\
        \n  -n, --numreads N\t\tnumber of files to extract\
        \n  -o, --outdir /DIR/TO/OUT\tdirectory of the output\
        \n  -s, --suffix EXT\t\tfile extension or suffix to extract (default=fast5)\
        \n  -t, --threads N\t\tnumber of threads to use (default=8)\n"
      exit
      shift 2
      ;;
    *) break
      ;;
  esac
done

# default suffix is "fast5"
if [ ! -n "$suffix" ];then
  echo "-s[--suffix] argument not given, using default "fast5""
  suffix=fast5
fi
if [ ! -n "$t" ];then
  t=8
fi

# check dir
if [ -e $outdir ]; then
  num=`find $outdir -name "*$suffix*" | wc -l`
  if [ "$num" -ge "$n" ];then
    echo "$outdir already has files, check the outdir"
    exit
  fi
else
  [ -e $outdir ]||mkdir $outdir
fi
outdir=`readlink -f $outdir`

cd $outdir
echo "output dir $PWD"

# execute code
echo "untarring from: $infile"
echo "$n files"
#make fofn file to feed into tar extract
fofn=$outdir/extract.fofn
tar -tf $infile | grep "$suffix" | head -n $n > $fofn
parallel --no-notice -j $t tar xvzf $infile --occurrence {} < $fofn
rm $fofn
echo "Finished with extracting $n files from $infile"
