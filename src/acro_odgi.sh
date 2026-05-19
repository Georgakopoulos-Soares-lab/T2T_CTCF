#!/bin/bash

# odgi sort \
#     -i ../res/misc/acro_short_arms.smooth.final.og \
#     -o ../res/misc/acro_short_arms.smooth.final.sorted.og \
#     -p Ygs \
#     -t 6 \
#     -O \
#     -P

# odgi layout \
#     -i ../res/misc/acro_short_arms.smooth.final.sorted.og \
#     -o ../res/misc/acro_short_arms.smooth.final.sorted.og.lay \
#     -T ../res/misc/acro.lay.tsv \
#     -t 6 \
#     -P

for ref in "chr13#0#p_arm" "chr14#0#p_arm" "chr15#0#p_arm" "chr21#0#p_arm" "chr22#0#p_arm"; do
    echo Processing: "$ref"
    chr=$(echo "$ref" | cut -d'#' -f1)
    
    odgi position \
        -i ../res/misc/acro_short_arms.smooth.final.sorted.og \
        -r "$ref" \
        -I \
        -b ../res/misc/acro_ctcf.bed \
        -t 6 \
        -P > "../res/misc/acro_short_arms.position.${chr}.tsv"

done
