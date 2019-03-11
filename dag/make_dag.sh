#!/bin/bash
dir=$(dirname "$0")
cd $dir/../
snakemake --dag parse_nanopore | dot -Tsvg > dag/nanopore_dag.svg
snakemake --dag nanopore_sv | dot -Tsvg > dag/nanopore_sv_dag.svg
snakemake --dag nanopore_methylation | dot -Tsvg > dag/nanopore_methylation_dag.svg
