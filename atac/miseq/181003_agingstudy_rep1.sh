#!/bin/bash
# As QC, I'm going to just go up to remove duplications
rawdir=/dilithium/Data/NGS/Raw/181003_chosigmaATACrep1_methylationStandards
root=/dilithium/Data/NGS/Aligned/181003_choSigmaATACrep1
btidx=/mithril/Data/NGS/Reference/cho/criGri1/criGri1

# get sample names
fqs=$(find $rawdir -name "cho*R1*fastq.gz")
for fq in $fqs;do
  base=$(basename "$fq")
  samp=${base%%_*}
  samps="$samps $samp"
done

trimdir=$root/fqtrim
[ -e $trimdir ]||mkdir -p $trimdir
if [ "$1" == "trim" ];then
  echo "trim"
  log=$trimdir/{}.trim.log
  com="trim_galore -v > $log &&\
    trim_galore --paired --fastqc \
    $rawdir/{}*R1*fastq.gz \
    $rawdir/{}*R2*fastq.gz \
    -o $trimdir &> $log"
  parallel "$com" ::: $samps
fi

bamdir=$root/bam
[ -e $bamdir ]||mkdir $bamdir
if [ "$1" == "align" ];then
  for samp in $samps;do
    echo $samp
    bam=$bamdir/$samp.sorted.bam
    log=$bamdir/$samp.align.log
    fq1=$(find $trimdir -name "$samp*R1*fq.gz")
    fq2=$(find $trimdir -name "$samp*R2*fq.gz")
    bowtie2 --version @> $log
    bowtie2 -k 4 -X 2000 -p 8 -t --local \
      -x $btidx -1 $fq1 -2 $fq2 \
      2> $log |\
      samtools sort -T $bamdir/$samp.sorted -o $bam -
    samtools index $bam
  done
fi

if [ "$1" == "bamprocess" ];then
  tmpdir=$root/tmp/bamprocess
  [ -e $tmpdir ]||mkdir -p $tmpdir
  bam=$bamdir/{}.sorted.bam
  multimap=$tmpdir/{}.multi.sam
  outbam=$bamdir/{}.processed.bam
#  parallel "samtools flagstat $bam &> $bamdir/{}.flagstst.txt" ::: $samps
  assigner="/home/isac/Code/projects/ambic-epigenome/atac/script/assign_multimappers.py"
  com="samtools view -F 524 -f 2 -u $bam |\
    samtools sort -@ 8 -T $tmpdir/{}.tmp -n -O sam - |\
    $assigner -k 4 --paired-end > $multimap &&\
    samtools fixmate -r $multimap -|\
    samtools view -F 1804 -f 2 -u - |\
    samtools sort -T $tmpdir/{}.sorted -@ 8 -o $outbam -"
  echo $com
#  parallel "$com" ::: $samps
  # picard
  jopt="-Xmx2G -Xms256M -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$tmpdir"
  for samp in $samps;do
    echo $samp
    bam=$bamdir/$samp.processed.bam
    dupmark=$bamdir/$samp.dupmark.bam
#    picard MarkDuplicates $jopt \
#      INPUT=$bam OUTPUT=$dupmark \
#      METRICS_FILE=$bamdir/$samp.duplicate_metrics.txt \
#      VALIDATION_STRINGENCY=LENIENT \
#      ASSUME_SORTED=true REMOVE_DUPLICATES=false \
#      &> $bamdir/$samp.markdup.log
#    bedtools bamtobed -i $dupmark | \
#      awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
#      grep -v 'chrM' | sort | uniq -c | \
#      awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} \
#    ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
#    END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; \
#    printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n" \
#    ,mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' \
#      > ${bamdir}/${samp}_pbcqc.txt
    nodup=$bamdir/$samp.nodup.bam
    samtools view -F 1804 -f 2 -b $dupmark > $nodup
    samtools index $nodup
    samtools flagstat $nodup > $bamdir/$samp.nodup.flagstat.txt
  done
fi

if [ "$1" == "qc" ];then
  qcdir=$root/qc
  qcout=$qcdir/readnumber_qc.txt
  echo "sample fastq bam bam-mito" | tr " " "\t" > $qcout
  for samp in $samps;do
    echo $samp
    fq=$(find $rawdir -name "${samp}_*R1*.fastq.gz")
    fqlines=$(gunzip -c $fq | wc -l)
    bam=$bamdir/$samp.nodup.bam
    bampairs=$(bedtools bamtobed -i $bam |\
      wc -l)
    nomito=$(bedtools bamtobed -i $bam |\
      grep -v "chrM" | wc -l)
    echo "$samp $fqlines $bampairs $nomito" |\
      awk '{ $2=$2/4;$3=$3/2;$4=$4/2; print $0,$4/$2 }' |\
      tr " " "\t" >> $qcout
  done
fi
