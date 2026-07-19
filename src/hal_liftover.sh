#!/bin/bash
hal_liftover() {
  local input_bed="$1"
  local output_bed="$2"

  halLiftover \
    --hdf5InMemory \
    ../data/multiple_sequence_alignments/8-t2t-apes-2023v2.hal \
    hs1 \
    "$input_bed" \
    hg38 \
    "$output_bed"
}

for file in chm13_MA0139.2.bed chm13_MA1929.2.bed chm13_MA1930.2.bed; do
  hal_liftover \
    "../res/misc/${file}" \
    "../res/misc/${file%.bed}.liftover.hg38.bed"
done