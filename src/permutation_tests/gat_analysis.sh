#!/bin/bash

GAT_dir="../../res/misc/gat"
log_dir="$GAT_dir/logs"

mkdir -p $log_dir

gat-run.py \
    --segment-file="$GAT_dir"/ctcf.bed \
    --workspace-file="$GAT_dir"/workspace.bed \
    --annotation-file="$GAT_dir"/censat.bed \
    --isochore-bed-file="$GAT_dir"/isochores.bed \
    --counter=segment-overlap \
    --num-samples=100 \
    --num-threads=4 \
    --random-seed=42 \
    --log="$log_dir"/gat.log \
    --with-segment-tracks \
    --verbose=5 > "$GAT_dir"/GAT_results.tsv