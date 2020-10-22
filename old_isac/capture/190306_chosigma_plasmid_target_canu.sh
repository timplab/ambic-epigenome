#!/bin/bash
dir=/kyber/Data/Nanopore/projects/ambic/sigma/reads/190225_choSigma_plasmid_target
for bam in $(find $dir -name "*bam"); do
  base=$(basename "$bam")
  base=${base%%.*}
  fq=$dir/${base}_sigmaIgG_insertion.fastq.gz
  if [ "$1" == "mkfq" ];then
  vcf=$dir/$base.sniffles.vcf
    qnames=$dir/$base.insert.qnames.txt
    samtools view -h $bam ambic_sigma_IgG_HC ambic_sigma_IgG_LC | cut -f1 > $qnames
    grep sigma $vcf | cut -d";" -f8 |\
      sed "s/RNAMES=//" |\
      tr "," "\n" | sort >> $qnames
    insertsam=$dir/$base.insert.sam
    sam=$dir/$base.sam
    samtools view -H $bam > $insertsam
    samtools view $bam > $sam 
    awk 'NR==FNR{c[$1]++;next};c[$1]' $qnames $sam >> $insertsam
    samtools fastq $insertsam | gzip > $fq
  fi
  if [ "$1" == "enrichment" ];then
    all=$(($(cat $dir/$base.fastq | wc -l)/4))
    target=$(($(gunzip -c $fq | wc -l)/4))
    echo $all,$target,$(($target*100/$all))%
  fi
done
for samp in CHOZNUnstableNoglnDay0 CHOZNStableNoglnDay0; do
  fq=$dir/${samp}_sigmaIgG_insertion.pooled.fastq.gz
  if [ "$1" == "merge" ];then
    fqs=$(find $dir -name "$samp*insertion.fastq.gz")
    gunzip -c $fqs | gzip > $fq
  fi
  if [ "$1" == "canu" ];then
    echo "canuing $samp"
    canu=/home/isac/Code/canu/Linux-amd64/bin/canu
    outdir=$dir/canu
    [ -e $outdir ]||mkdir $outdir
    log=$outdir/$samp.canu.log
    $canu \
      -p $samp -d $outdir/$samp \
      -genomeSize=10k \
      -nanopore-raw $fq &> $log
  fi
  pre=$dir/canu/$samp/$samp
  if [ "$1" == "align" ];then
    echo "aligning $samp"
    idx=/mithril/Data/NGS/Reference/cho/chok1/genome/Cricetulus_griseus_chok1gshd.sigmaIgG.mmi
    pre=$dir/canu/$samp/$samp
    fa=$pre.contigs.fasta
    out=$pre.sorted.bam
    log=$pre.align.log
    minimap2 --MD -t 10 -ax map-ont -L $idx $fa 2> $log |\
      samtools view -b -q 20 - |\
      samtools sort -o $out
    samtools index $out
  fi
  if [ "$1" == "locate" ];then
    echo $samp
    out=$pre.contig.summary.txt
    python ./insert_locater.py -b $pre.sorted.bam \
      -c ambic_sigma_IgG_HC ambic_sigma_IgG_LC > $out
  fi
done
