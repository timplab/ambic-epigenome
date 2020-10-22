#!/bin/bash
export LC_ALL=C
binroot="/mnt/ambic1"
export PATH="$binroot/bin:$binroot/Code/miniconda3/bin:$PATH"
root=/shared/Data
calldir=$root/bcall
bamdir=$root/bam
mroot=$root/mcall
mbeddir=$root/mbed

# for bcall, use bash scripts
if [ "$1" == "bcall" ];then
  echo "bcall"
  logdir=
fi

# snakemake
snake=ambic_dataparse_snakemake.py
config=ambic_snakemake.config
clustercfg=ambic_cluster_config.json
if [ "$1" == "align" ];then
  logdir=/shared/Data/log/guppy_basecall
  [ -e $logdir ]||mkdir -p $log
  snakemake -j 25 \
    --snakefile $snake --configfile $config \
    --cluster-config $clustercfg \
    --cluster "sbatch -o {cluster.log} -e {cluster.log} \
    -c {cluster.cores} -J {cluster.jobname}" \
    $bams -p
fi
