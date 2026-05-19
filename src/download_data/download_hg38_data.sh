#!/bin/bash
set -euo pipefail

# GRCh38.p14 downloads
mkdir -p ../../data/hg38
cd ../../data/hg38

# Create an array with the URLs of the metadata needed from GRCh38.p14
hg38_urls=(
    # Genome FASTA
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    
    # GFF3 files with gene annotations
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
    
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt 
)

printf '%s\n' "${hg38_urls[@]}" | parallel -j 2 wget

echo "Decompressing .gz files"
find . -name "*.gz" | parallel -j 2 gunzip {}

echo "Compressing and indexing FASTA"

bgzip -@ 8 GCF_000001405.40_GRCh38.p14_genomic.fna
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna.gz
cut -f1,2 GCF_000001405.40_GRCh38.p14_genomic.fna.gz.fai > GCF_000001405.40_GRCh38.p14_genomic.chrom.sizes