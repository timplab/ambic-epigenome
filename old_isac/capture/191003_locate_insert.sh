#!/bin/bash
# chozn
dir=/home/isac/Data/ambic/targeted/190225_chosigma_targeted
bam=$(find $dir -name "*sorted.bam")
echo $bam
base=$(basename "$bam")
base=${base%%.*}
#python insert_locater.py -b $bam -c ambic_sigma_IgG_HC ambic_sigma_IgG_LC
# nih
dir=/home/isac/Data/ambic/targeted/191001_CHO_NIH_IgG3_targeted
bam=$(find $dir -name "*sorted.bam")
echo $bam
base=$(basename "$bam")
base=${base%%.*}
python insert_locater.py -b $bam -c 4_0cdhfr_vrc01wtg1m3_dgv 
