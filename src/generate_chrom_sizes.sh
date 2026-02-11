#!/bin/bash

find ../data/ -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fa" \) | \
  parallel -j 60 --bar 'samtools faidx {} && cut -f1,2 {}.fai > {}.chrom.sizes'