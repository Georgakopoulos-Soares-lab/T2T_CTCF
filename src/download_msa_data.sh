#!/bin/bash
# Multiple sequences alignment downloads for CHM13 vs T2T primates and CHM13 vs human pangenome assemblies
mkdir -p ../data/multiple_sequence_alignments
cd ../data/multiple_sequence_alignments

# Create an array with the URLs of the metadata needed from HG002
msa_urls=(
    # 8-way CACTUS alignment in HAL format for CHM13, hg38, and the 6 T2T primate species
    "https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/8-t2t-apes-2023v2/8-t2t-apes-2023v2.hal"

    # Phylogenetic tree for these species in NEWICK format
    "https://cgl.gi.ucsc.edu/data/cactus/t2t-apes/8-t2t-apes-2023v2/8-t2t-apes-2023v2.nh"

    # CACTUS alignment in HAL format for the version 1 human pangenome assemblies in HAL format 
    "https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.hal"
)

printf '%s\n' "${msa_urls[@]}" | parallel -j 4 wget