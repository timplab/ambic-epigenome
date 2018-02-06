#!/bin/bash
bamdir=/atium/Data/Nanopore/Analysis/170608_cho/hitbam_nometh/oldgenome
bampath=${bamdir}/170608_choIgGNIH_hits.sorted.bam
covpath=${bamdir}/170608_choIgGNIH_hits_coverage.txt

samtools depth $bampath > $covpath
