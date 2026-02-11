#!/bin/bash

# Parse command-line arguments
motif_file=$1 # CTCF .meme file
background_file=$2 # 1-based Markov background file
fasta_file=$3 # Genome .fasta file 
output_base=$4 # Output name

# Extract filenames for naming output files
filename=$(basename "${fasta_file}")
motif_name=$(basename ${motif_file} .meme)

# Create final output directory if it doesn't exist
final_output_dir="${output_base}"
mkdir -p ${final_output_dir}

# Create temporary directory for FIMO output
temp_output_dir="${output_base}/temp_fimo_res_${motif_name}_$$"
mkdir -p ${temp_output_dir}

# Run FIMO motif scanning
fimo \
    --max-stored-scores 1000000000 \
    --thresh 1e-6 \
    --bfile "$background_file" \
    --oc ${temp_output_dir} ${motif_file} ${fasta_file}

# Move the main results file to final location with descriptive name
mv "${temp_output_dir}/fimo.tsv" "${final_output_dir}/${motif_name}.tsv"

# Clean up temporary directory and intermediate files
rm -rf ${temp_output_dir}