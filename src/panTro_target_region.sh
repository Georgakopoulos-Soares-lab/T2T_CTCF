#!/bin/bash
set -euo pipefail

DATA=../data
OUT=../res/misc
THREADS=80
mkdir -p "$OUT"

GENOMES=(human pan_pan gor_gor pon_abe pon_pyg sym_syn)

declare -A REF
REF[human]="$DATA/chm13/chm13v2.0.fa.gz"
REF[pan_pan]="$DATA/primates/PanPan/GCF_029289425.2_NHGRI_mPanPan1-v2.1_pri_genomic.fna.gz"
REF[gor_gor]="$DATA/primates/GorGor/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna.gz"
REF[pon_abe]="$DATA/primates/PonAbe/GCF_028885655.2_NHGRI_mPonAbe1-v2.1_pri_genomic.fna.gz"
REF[pon_pyg]="$DATA/primates/PonPyg/GCF_028885625.2_NHGRI_mPonPyg2-v2.1_pri_genomic.fna.gz"
REF[sym_syn]="$DATA/primates/SymSyn/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna.gz"

QUERY_SRC="$DATA/primates/PanTro/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz"

# QUERY_REGION="chr1_hap1_hsa1:209537154-209792028"
QUERY_REGION="chr1_hap1_hsa1:1-231299448"

samtools faidx "$QUERY_SRC" "$QUERY_REGION" > "$OUT/panTro_target_region.fa"

sed "s/^>/>panTro#/" "$OUT/panTro_target_region.fa" > "$OUT/panTro_target_region.tagged.fa"
QUERY="$OUT/panTro_target_region.tagged.fa"

PAF="$OUT/regions_vs_primates.paf"
: > "$PAF"

for g in "${GENOMES[@]}"; do
  echo ">> mapping query vs $g" >&2
  minimap2 -x asm20 -c --eqx -p 0.5 -N 50 -t "$THREADS" \
    <(zcat "${REF[$g]}" | sed "s/^>/>${g}#/") \
    "$QUERY" >> "$PAF"
done