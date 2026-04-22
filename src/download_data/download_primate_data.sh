#!/bin/bash
set -euo pipefail

# T2T primate downloads
mkdir -p ../../data/primates
cd ../../data/primates

mkdir -p SymSyn PonPyg PonAbe GorGor PanPan PanTro

# Create an array with the URLs of the metadata needed from T2T Primates project
primate_urls=(
    # Symphalangus syndactylus data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.1_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.0_pri.RepeatMasker_v1.1.bb" # Maybe switch chr12 and chr19??
    "https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/epi/mSymSyn1_v2.1.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Pongo pygmaeus data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_pri.RepeatMasker_v1.1.bb"
    "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/epi/mPonPyg2_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Pongo abelii data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_pri.RepeatMasker_v1.1.bb"
    "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/epi/mPonAbe1_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Gorilla gorilla data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_pri.RepeatMasker_v1.1.bb"
    "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/epi/mGorGor1_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Pan paniscus data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_pri.RepeatMasker_v1.1.bb"
    "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/epi/mPanPan1_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Pan troglodytes data
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gff.gz"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_assembly_report.txt"
    "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_pri.RepeatMasker_v1.1.bb"
    "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v2.0.bb"
    "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/epi/mPanTro3_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw"

    # Conservation scores of homo sapiens compared to the primate assemblies
    "https://garrisonlab.s3.amazonaws.com/t2t-primates/wfmash-v0.13.0/conservation/approach2_with_cds/chm13_1/scores/approach2.scores_v0.3.bw"
)

printf '%s\n' "${primate_urls[@]}" | parallel -j 12 wget

pigz -d *.gz

for bb in *.bb; do
    bigBedToBed "$bb" "${bb%.bb}.bed"
done

for fa in *.fna; do
    base="${fa%_genomic.fna}"
    report="${base}_assembly_report.txt"
    
    seqkit replace -p '^(\S+)' -r '{kv}' \
        -k <(grep -v '^#' "$report" | awk -F'\t' '$7 != "na" {print $7"\t"$1}') \
        "$fa" -o "${fa}.tmp" \
    && mv "${fa}.tmp" "$fa"
done

for gff in *.gff; do
    base="${gff%_genomic.gff}"
    report="${base}_assembly_report.txt"

    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==FNR {
            if ($0 !~ /^#/ && $7 != "na") map[$7] = $1
            next
        }

        /^##sequence-region/ {
            n = split($0, a, /[ \t]+/)
            if (a[2] in map) a[2] = map[a[2]]
            line = a[1]
            for (i=2; i<=n; i++) line = line " " a[i]
            print line
            next
        }

        /^#/ { print; next }

        {
            if ($1 in map) $1 = map[$1]
            print
        }
    ' "$report" "$gff" > "${gff}.tmp" && mv "${gff}.tmp" "$gff"

    gff3sort.pl --precise "$gff" > "${base}.sorted.gff"

    mv "${base}.sorted.gff" "$gff"

    bgzip -@ 8 "$gff"
    tabix -p gff "${gff}.gz"
done

for fa in *.fna; do
  base="${fa%.fna}"
  bgzip -@ 8 "$fa" && \
  samtools faidx "${base}.fna.gz" && \
  cut -f1,2 "${base}.fna.gz.fai" > "${base}.chrom.sizes"
done

species=(
    "GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri|mSymSyn1|SymSyn"
    "GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri|mPonPyg2|PonPyg"
    "GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri|mPonAbe1|PonAbe"
    "GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri|mGorGor1|GorGor"
    "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri|mPanPan1|PanPan"
    "GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri|mPanTro3|PanTro"
)

for entry in "${species[@]}"; do
    IFS='|' read -r ncbi ga dir <<< "$entry"
    mv \
        "${ncbi}_genomic.fna.gz" \
        "${ncbi}_genomic.fna.gz.fai" \
        "${ncbi}_genomic.fna.gz.gzi" \
        "${ncbi}_genomic.chrom.sizes" \
        "${ncbi}_genomic.gff" \
        "${ncbi}_assembly_report.txt" \
        ${ga}_*.bed \
        ${ga}_*.bb \
        ${ga}_*.bw \
        "$dir/"
done
