#!/bin/bash
# T2T primate downloads
mkdir -p ../data/primates
cd ../data/primates

mkdir -p SymSyn PonPyg PonAbe GorGor PanPan PanTro

# Create an array with the URLs of the metadata needed from T2T Primates project
primate_urls=(
    # # Symphalangus syndactylus data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/878/055/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri/GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.1_CenSat_v2.0.bb"
    # "https://genomeark.s3.amazonaws.com/species/Symphalangus_syndactylus/mSymSyn1/assembly_curated/repeats/mSymSyn1_v2.0_pri.RepeatMasker_v1.1.bb" # Maybe switch chr12 and chr19??

    # # Pongo pygmaeus data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/625/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_CenSat_v2.0.bb"
    # "https://genomeark.s3.amazonaws.com/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/repeats/mPonPyg2_v2.0_pri.RepeatMasker_v1.1.bb"
    
    # # Pongo abelii data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/885/655/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri/GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_CenSat_v2.0.bb"
    # "https://genomeark.s3.amazonaws.com/species/Pongo_abelii/mPonAbe1/assembly_curated/repeats/mPonAbe1_v2.0_pri.RepeatMasker_v1.1.bb"

    # # Gorilla gorilla data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_pri.RepeatMasker_v1.1.bb"
    # "https://genomeark.s3.amazonaws.com/species/Gorilla_gorilla/mGorGor1/assembly_curated/repeats/mGorGor1_v2.0_CenSat_v2.0.bb"

    # # Pan paniscus data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_CenSat_v2.0.bb"
    # "https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0_pri.RepeatMasker_v1.1.bb"

    # # Pan troglodytes data
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gff.gz"
    # "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_assembly_report.txt"
    # "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_pri.RepeatMasker_v1.1.bb"
    # "https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0_CenSat_v2.0.bb"

    # Conservation scores of homo sapiens compared to the primate assemblies
    https://garrisonlab.s3.amazonaws.com/t2t-primates/wfmash-v0.13.0/conservation/approach2_with_cds/chm13_1/scores/approach2.scores_v0.3.bw
)

printf '%s\n' "${primate_urls[@]}" | parallel -j 12 wget

mv approach2.scores_v0.3.bw conservation_scores.bw

# pigz -d *.gz

# for bb in *.bb; do
#     bigBedToBed "$bb" "${bb%.bb}.bed"
# done

# mv \
#     GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.fna \
#     GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_genomic.gff \
#     GCF_028878055.3_NHGRI_mSymSyn1-v2.1_pri_assembly_report.txt \
#     mSymSyn1_*.bed \
#     mSymSyn1_*.bb \
#     SymSyn/

# mv \
#     GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna \
#     GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.gff \
#     GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_assembly_report.txt \
#     mPonPyg2_*.bed \
#     mPonPyg2_*.bb \
#     PonPyg/

# mv \
#     GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.fna \
#     GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_genomic.gff \
#     GCF_028885655.2_NHGRI_mPonAbe1-v2.0_pri_assembly_report.txt \
#     mPonAbe1_*.bed \
#     mPonAbe1_*.bb \
#     PonAbe/

# mv \
#     GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.fna \
#     GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_genomic.gff \
#     GCF_029281585.2_NHGRI_mGorGor1-v2.1_pri_assembly_report.txt \
#     mGorGor1_*.bed \
#     mGorGor1_*.bb \
#     GorGor/

# mv \
#     GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna \
#     GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gff \
#     GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_assembly_report.txt \
#     mPanPan1_*.bed \
#     mPanPan1_*.bb \
#     PanPan/

# mv \
#     GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna \
#     GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.gff \
#     GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_assembly_report.txt \
#     mPanTro3_*.bed \
#     mPanTro3_*.bb \
#     PanTro/