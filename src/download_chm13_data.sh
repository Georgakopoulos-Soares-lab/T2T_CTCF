#!/bin/bash
# CHM13 downloads
mkdir -p ../data/chm13
cd ../data/chm13

# Create an array with the URLs of the metadata needed from CHM13
chm13_urls=(
    # Genome FASTA
    # CHM13 v2.0 assembly (AC:GCA_009914755.4) with sequences soft-masked by the T2T team and sequence names converted to chr*. 
    # chrY in this file was assembled from sample HG002/NA24385. The original sequence accession numbers are shown in the FASTA header.
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
    
    # GFF3 file with gene annotations
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz"
    
    # Cytobands
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed"

    # Segmental duplications in simple format
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed"
    
    # CenSat annotations 
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed"
    
    # Repeat masker annotations
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
    
    # Composite repeats
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_composite-repeats_2022DEC.bed"

    # New satellites
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_new-satellites_2022DEC.bed"
    
    # X/Y sequence class
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_chrXY_sequence_class_v1.bed"
    
    # Telomere annotations
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed"

    # Y specific annotations
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_inverted_repeats_v1.bed"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_amplicons_v1.bed"
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_AZF_DYZ_v1.bed"

    # CpG methylation levels of CHM13 
    "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/regulation/chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw"

)

printf '%s\n' "${chm13_urls[@]}" | parallel -j 8 wget

pigz -d chm13v2.0.fa.gz # Decompress FASTA file
pigz -d chm13v2.0_RefSeq_Liftoff_v5.2.gff3.gz # Decompress GFF3 file