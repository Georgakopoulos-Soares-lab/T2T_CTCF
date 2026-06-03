#!/bin/bash

GAT_dir="../../res/misc/gat"
log_dir="$GAT_dir/logs"
per_chrom_results_dir="$GAT_dir/per_chrom_results"
primate_results_dir="$GAT_dir/primates"

mkdir -p "$log_dir" "$per_chrom_results_dir" "$primate_results_dir"

# # Genome-wide analysis
# gat-run.py \
#     --segment-file="$GAT_dir"/ctcf.bed \
#     --workspace-file="$GAT_dir"/workspace.bed \
#     --annotation-file="$GAT_dir"/censat.bed \
#     --isochore-bed-file="$GAT_dir"/isochores.bed \
#     --counter=segment-overlap \
#     --num-samples=10000 \
#     --num-threads=140 \
#     --random-seed=42 \
#     --log="$log_dir"/gat.log \
#     --with-segment-tracks \
#     --verbose=5 > "$GAT_dir"/GAT_results.tsv

# # Analysis per chromosome
# chroms=$(cut -f1 "$GAT_dir/workspace.bed" | sort -u)

# for chrom in $chroms; do
#     gat-run.py \
#         --segment-file=<(awk -v c="$chrom" '$1==c' "$GAT_dir/ctcf.bed") \
#         --workspace-file=<(awk -v c="$chrom" '$1==c' "$GAT_dir/workspace.bed") \
#         --annotation-file=<(awk -v c="$chrom" '$1==c' "$GAT_dir/censat.bed") \
#         --isochore-bed-file=<(awk -v c="$chrom" '$1==c' "$GAT_dir/isochores.bed") \
#         --counter=segment-overlap \
#         --num-samples=10000 \
#         --num-threads=140 \
#         --random-seed=42 \
#         --log="$log_dir/gat_${chrom}.log" \
#         --with-segment-tracks \
#         --verbose=5 > "$per_chrom_results_dir/GAT_results_${chrom}.tsv"
# done

# # Run analysis for gamma-satellites
# gat-run.py \
#     --segment-file="$GAT_dir"/ctcf.bed \
#     --workspace-file="$GAT_dir"/workspace.bed \
#     --annotation-file="$GAT_dir"/gsat.bed \
#     --isochore-bed-file="$GAT_dir"/isochores.bed \
#     --counter=segment-overlap \
#     --num-samples=10000 \
#     --num-threads=140 \
#     --random-seed=42 \
#     --log="$log_dir"/gat_gsat.log \
#     --with-segment-tracks \
#     --verbose=5 > "$GAT_dir"/GAT_results_gsat.tsv

# # Run analysis for pseudo-homologous regions
# gat-run.py \
#     --segment-file="$GAT_dir"/phr_analysis/ctcf_sites.bed \
#     --workspace-file="$GAT_dir"/phr_analysis/workspace.bed \
#     --annotation-file="$GAT_dir"/phr_analysis/phrs.bed \
#     --isochore-bed-file="$GAT_dir"/phr_analysis/isochores.bed \
#     --counter=segment-overlap \
#     --num-samples=10000 \
#     --num-threads=6 \
#     --random-seed=42 \
#     --with-segment-tracks \
#     --log="$log_dir"/gat.log \
#     --verbose=5 > "$GAT_dir"/phr_analysis/GAT_results.tsv

# Run analysis for primate species
for species in chm13 gorgor panpan pantro ponabe ponpyg symsyn; do
    
    echo "Running GAT for $species..."
    
    gat-run.py \
        --segment-file="$primate_results_dir"/"$species"/ctcf_sites.bed \
        --workspace-file="$primate_results_dir"/"$species"/workspace.bed \
        --annotation-file="$primate_results_dir"/"$species"/censat.bed \
        --isochore-bed-file="$primate_results_dir"/"$species"/isochores.bed \
        --counter=segment-overlap \
        --num-samples=10 \
        --num-threads=5 \
        --random-seed=42 \
        --log="$log_dir"/"$species"_gat.log \
        --with-segment-tracks \
        --verbose=5 > "$primate_results_dir"/"$species"/GAT_results.tsv
done
