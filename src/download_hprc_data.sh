#!/bin/bash
# Human Pangenome Resource download
mkdir -p ../data/human_pangenome
cd ../data/human_pangenome

awk -F'\t' '{print $1}' ../../links.txt | sort -u | xargs -I {} mkdir -p {}

parallel -j 16 --colsep '\t' 'aws s3 --no-sign-request cp {2} {1}/' :::: ../../links.txt

echo "Decompressing .gz files"
find . -name "*.gz" | parallel -j 64 gunzip {}