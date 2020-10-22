#!/bin/bash
qcdir=/mnt/c/qc
date=190320
qc=$qcdir/${date}_fastq_qc.csv
sum=$qcdir/${date}_ambic_agingstudy_fastq_summary.csv
scr=./ambic_agingstudy_fastq_qc.py

python $scr $qc > $sum
