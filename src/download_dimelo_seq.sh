#!/bin/bash
# CHM13 downloads
mkdir -p ../data/dimelo_seq
cd ../data/dimelo_seq

# Create an array with the URLs of the metadata needed from CHM13
dimelo_seq_urls=(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6337nnn/GSM6337767/suppl/GSM6337767%5Fgm12878%5Fh3k27ac%2Emod%5Fmappings%2Esorted%2Ebam"
    
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6337nnn/GSM6337768/suppl/GSM6337768%5Fgm12878%5Fh3k27me3%2Emod%5Fmappings%2Esorted%2Ebam"
    
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6337nnn/GSM6337769/suppl/GSM6337769%5Fgm12878%5Fh3k4me3%2Emod%5Fmappings%2Esorted%2Ebam"
)

printf '%s\n' "${dimelo_seq_urls[@]}" | parallel -j 8 wget
