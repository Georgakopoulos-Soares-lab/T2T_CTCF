#!/bin/bash
# HG002 downloads
mkdir -p ../data/hg002
cd ../data/hg002

# Create an array with the URLs of the metadata needed from HG002
hg002_urls=(
    # Genome FASTA paternal and maternal
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.pat.fasta.gz"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.mat.fasta.gz"
    
    # GFF3 files with gene annotations for paternal and maternal
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/JHULiftoff/v0.6/hg002.v1.1.loff.v0.6.mat.gff.gz"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/JHULiftoff/v0.6/hg002.v1.1.loff.v0.6.pat.gff.gz"
    
    # Cytobands
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/cytoband/cytoBand.hg002v1.1.mat.lifted_from_v1.0.1.bb"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/cytoband/cytoBand.hg002v1.1.pat.lifted_from_v1.0.1.bb"

    # Segmental duplications in simple format
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/segdups/hg002v1.1.pat.SDs.013025.bed"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/segdups/hg002v1.1.mat.SDs.013025.bed"
    
    # CenSat annotations 
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.MAT.bed"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/centromere/hg002v1.1_v2.0/hg002v1.1.cenSatv2.0.PAT.bed"
    
    # Repeat masker annotations
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/repeatmasker/hg002v1.1.maternal.bb"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/repeatmasker/hg002v1.1.paternal.bb"

    # CpG methylation levels of HG002
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/annotation/methylation/HG002v1.1.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # GpC methylation in Z-score format from NanoNOMe data for the HG002 cell line 
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/regulation/chm13v2.0_hg002_GpC_zscore_ont_guppy5.0.7_nanopolish0.13.2.bw"
)

printf '%s\n' "${hg002_urls[@]}" | parallel -j 8 wget

echo "Decompressing .gz files"
find . -name "*.gz" | parallel -j 64 gunzip {}

for bb in *.bb; do
    bigBedToBed "$bb" "${bb%.bb}.bed"
done