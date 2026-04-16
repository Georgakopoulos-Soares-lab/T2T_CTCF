#!/bin/bash
#SBATCH -J modkit_prepare
#SBATCH -o ../logs/modkit_prepare.out
#SBATCH -e ../logs/modkit_prepare.err
#SBATCH -p gg
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -A MCB25091

source ~/.bashrc
conda activate dimelo-toolkit

cd ../data/dimelo_seq

for bam in *.sorted.bam; do
    out="${bam%.sorted.bam}.updated.bam"
    echo "Processing $bam -> $out"
    modkit update-tags \
        --mode ambiguous \
        --threads 130 \
        "$bam" \
        "$out"

    samtools index --threads 130 "$out"
done

cd ../chm13

awk -v OFS='\t' '$1 !~ /^(chrM|chrY)$/ {print $1, 0, $2}' chm13v2.0.fa.fai > chm13v2.0.genome.bed