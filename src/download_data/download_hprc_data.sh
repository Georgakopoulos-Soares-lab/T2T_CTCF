#!/bin/bash
set -euo pipefail

# Human Pangenome Resource download
mkdir -p ../../data/human_pangenome
cd ../../data/human_pangenome

awk -F'\t' '{print $1}' ../../links.txt | sort -u | xargs -I {} mkdir -p {}

parallel -j 16 --colsep '\t' 'aws s3 --no-sign-request cp {2} {1}/' :::: ../../links.txt

echo "Decompressing .gz files"
find . -name "*.gz" | parallel -j 64 gunzip {}

wget "https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-chm13.wave.vcf.gz"
wget "https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-chm13.wave.vcf.gz.tbi"
wget "https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-chm13.gfa.gz"

echo "Compressing and indexing FASTA files"

index_fasta() {
    local fa="$1"
    local base="${fa%.fa}"

    bgzip -@ 8 "$fa" && \
    samtools faidx "${base}.fa.gz" && \
    cut -f1,2 "${base}.fa.gz.fai" > "${base}.chrom.sizes"
}

export -f index_fasta

find . -type f -name "*.fa" -print0 \
  | parallel -0 -j 8 index_fasta {}