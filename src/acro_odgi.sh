#!/bin/bash

data_dir="../data/chm13"
res_dir="../res/misc"

bedtools getfasta -fi $data_dir/chm13v2.0.fa -bed $res_dir/acrocentric_arms.bed | \
  awk '/^>/{split($0,a,":"); chr=substr(a[1],2); print ">"chr"#0#p_arm"; next}{print}' \
  > $res_dir/acrocentric_arms.bedacro_short_arms.fa

bgzip $res_dir/acrocentric_arms.bedacro_short_arms.fa
samtools faidx $res_dir/acrocentric_arms.bedacro_short_arms.fa.gz
mkdir -p $res_dir/pggb_out

pggb \
    --input-fasta $res_dir/acro_short_arms.fa.gz \
    --segment-length 25000 \
    --map-pct-id 95 \
    --n-haplotypes 5 \
    1 100 \
    --threads 120 \
    --output-dir $res_dir/pggb_out \
    --multiqc \
    --poa-threads 40 \
    --compress

odgi sort \
    -i $res_dir/acro_short_arms.smooth.final.og \
    -o $res_dir/acro_short_arms.smooth.final.sorted.og \
    -p Ygs \
    -t 6 \
    -O \
    -P

odgi layout \
    -i $res_dir/acro_short_arms.smooth.final.sorted.og \
    -o $res_dir/acro_short_arms.smooth.final.sorted.og.lay \
    -T $res_dir/acro.lay.tsv \
    -t 6 \
    -P

for ref in "chr13#0#p_arm" "chr14#0#p_arm" "chr15#0#p_arm" "chr21#0#p_arm" "chr22#0#p_arm"; do
    echo Processing: "$ref"
    chr=$(echo "$ref" | cut -d'#' -f1)
    
    odgi position \
        -i $res_dir/acro_short_arms.smooth.final.sorted.og \
        -r "$ref" \
        -I \
        -b $res_dir/acro_ctcf.bed \
        -t 6 \
        -P > "$res_dir/acro_short_arms.position.${chr}.tsv"

done
