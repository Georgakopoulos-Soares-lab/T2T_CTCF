# Analysis of the potential CTCF binding sites across the complete telomere-to-telomere human genome

Comparative analysis of CTCF transcription factor binding sites across Telomere-to-Telomere (T2T) genome assemblies, including human references (CHM13, hg38), the HG002 diploid genome, 40+ human pangenome assemblies, and six T2T great ape genomes.

## Overview

[CTCF](https://en.wikipedia.org/wiki/CTCF) (CCCTC-binding factor) is a key transcription factor involved in chromatin organization and gene regulation. This pipeline performs genome-wide *in silico* identification of CTCF binding motifs using [FIMO](https://meme-suite.org/meme/doc/fimo.html) (Find Individual Motif Occurrences) and compares the resulting binding site landscapes across multiple genome assemblies to identify conserved and lineage-specific sites.

Three JASPAR CTCF motif variants are scanned:
- **MA0139.2** -- canonical CTCF motif
- **MA1929.2** -- CTCF variant motif
- **MA1930.2** -- CTCF variant motif

## Genomes analyzed

| Category | Assembly | Description |
|---|---|---|
| **Human references** | T2T-CHM13 v2.0 | Complete T2T human reference genome |
| | GRCh38.p14 | Standard human reference (hg38) |
| **Great apes** | *Pan troglodytes* | Chimpanzee (mPanTro3 v2.0) |
| | *Pan paniscus* | Bonobo (mPanPan1 v2.0) |
| | *Gorilla gorilla* | Western gorilla (mGorGor1 v2.1) |
| | *Pongo abelii* | Bornean orangutan (mPonAbe1 v2.0) |
| | *Pongo pygmaeus* | Sumatran orangutan (mPonPyg2 v2.0) |
| | *Symphalangus syndactylus* | Siamang (mSymSyn1 v2.1) |

## Repository structure

```
T2T_CTCF/
├── config.yml                     # Path configuration (edit before running)
├── environment.yml                # Conda environment specification
├── src/
│   ├── download_data/             # Genome/annotation/data download scripts
│   │   ├── download_chm13_data.sh
│   │   ├── download_hg38_data.sh
│   │   ├── download_hg002_data.sh
│   │   ├── download_hprc_data.sh
│   │   ├── download_primate_data.sh
│   │   ├── download_msa_data.sh        # HAL multiple sequence alignments
│   │   └── download_dimelo_seq.sh      # DiMeLo-seq methylation data
│   ├── FIMO.smk                   # Snakemake workflow for FIMO motif scanning
│   ├── FIMO.sh                    # FIMO runner script (called by Snakemake)
│   ├── hal_liftover.sh            # Liftover CTCF coords between assemblies via HAL
│   ├── format_gff3_files.sh       # Standardize/annotate CHM13 gene GFF3
│   ├── acro_odgi.sh               # Build acrocentric short-arm pangenome graph (odgi)
│   ├── panTro_target_region.sh    # Extract chimpanzee target regions
│   ├── modkit_prepare.sh          # Prepare ONT modkit methylation calls
│   ├── dimelo_supp.R              # DiMeLo-seq supplementary figure
│   ├── figure1/                   # Scripts generating Figure 1 panels
│   │   ├── create_hg38_chm13_beds.R
│   │   ├── venn_diagram.R
│   │   ├── karyoplot.R
│   │   ├── density_plots.R
│   │   ├── gsat_density_and_methylation.R
│   │   ├── new_genes_plots.R
│   │   └── ajax_teucer_analysis.R
│   ├── figure2/                   # Scripts generating Figure 2 panels
│   │   ├── acrocentric_short_arms.R
│   │   ├── acrocentric_ctcf_sharing.R
│   │   ├── ctcf_sdi_correlation.R
│   │   ├── ctcf_phr_coverage.R
│   │   └── plot_rdna.R
│   ├── figure3/                   # Scripts generating Figure 3 panels
│   │   ├── paf_visualization.R
│   │   ├── censat_conserv.R
│   │   ├── censat_conserv_per_chrom.R
│   │   ├── censat_density_plots.R
│   │   ├── conservation_beeswarm.R
│   │   ├── lineage_specific_SDs.R
│   │   ├── spacer_te_analysis.R
│   │   └── subterm_analysis.R
│   ├── permutation_tests/         # GAT-based permutation/enrichment testing
│       ├── prepare_permutation_data.R
│       ├── prepare_gsat_permutation_data.R
│       ├── prepare_phr_permutation_data.R
│       ├── prepare_primate_data.R
│       ├── prepare_teucer_permutation_data.R
│       └── gat_analysis.sh
├── data/                     # Input data (not tracked in git)
└── res/                      # Results (not tracked in git)
```

## Pipeline

### 1. Data download

Download genome assemblies and annotations from their respective repositories:

```bash
cd src/download_data
bash download_chm13_data.sh
bash download_hg38_data.sh
bash download_hg002_data.sh
bash download_hprc_data.sh
bash download_primate_data.sh
bash download_msa_data.sh
bash download_dimelo_seq.sh
```

**Data sources:**
- CHM13: [T2T Consortium](https://github.com/marbl/CHM13) (AWS S3)
- hg38: [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Great apes: [NCBI](https://www.ncbi.nlm.nih.gov/) & [GenomeArk](https://genomeark.github.io/)
- Multiple sequence alignments: 8-way T2T ape CACTUS alignment in HAL format

### 2. CTCF motif scanning (FIMO)

The Snakemake workflow auto-detects all FASTA files under `data/` and all motif files under `data/motifs/`, then:

1. Creates 1st-order Markov background models per genome (`fasta-get-markov`)
2. Scans each genome for each CTCF motif using FIMO (p-value threshold: 1e-6)

```bash
cd src
snakemake -s FIMO.smk --cores <N>
```

Results are written as TSV files to `res/fimo_res/<genome>/<motif>.tsv`.

### 3. BED conversion and coordinate liftover

Convert FIMO results for CHM13 and hg38 to BED format and lift coordinates between assemblies using the HAL alignment:

```bash
Rscript figure1/create_hg38_chm13_beds.R
bash hal_liftover.sh
```

The liftover uses `halLiftover` to map CTCF site coordinates between assemblies (e.g. CHM13/hs1 to hg38, or to primate genomes) via the T2T CACTUS multiple sequence alignment.

### 4. Comparative analysis and visualization

Figure-generating scripts are grouped by figure under `src/figure1/`, `src/figure2/`, and `src/figure3/`, with `src/permutation_tests/` providing GAT-based enrichment/permutation testing that feeds several of these figures. Each script reads its inputs from `res/` (as produced by the steps above) and writes its output to `res/plots/`. Run them with `Rscript <path>` from `src/`, e.g. `Rscript figure1/karyoplot.R`.

#### Figure 1

**Karyotype plot** -- Visualize CTCF binding site density (per 1 Mb window) across the CHM13 karyotype, overlaid with cytobands.

```bash
Rscript figure1/karyoplot.R          # -> res/plots/karyotype_plot.pdf
```

**Venn diagrams** -- Compare CTCF sites between hg38 and T2T-CHM13 (via liftover). Sites are classified as shared (overlapping within the lifted coordinates) or assembly-specific. Regions that fail liftover are counted as CHM13-unique.

```bash
Rscript figure1/venn_diagram.R       # -> res/plots/all_liftover_venns.pdf
```

**Density plots** -- Calculate CTCF density (number of sites per 1Mb) across different CenSat types and chromosomes. The plots created depict CTCF density in different CenSat compartments compared to genome wide for all 3 JASPAR CTCF motifs used.

```bash
Rscript figure1/density_plots.R       # -> res/plots/combined_densities.svg
```

**Gamma-satellite density and methylation** -- Barplot of CTCF density across gamma-satellite (GSAT) families and a boxplot of the distribution of methylation scores across those families.

```bash
Rscript figure1/gsat_density_and_methylation.R   # -> res/plots/gsat_plus_meth.svg
```

**New genes analysis** -- Summarizes CTCF sites overlapping newly annotated CHM13 genes (stacked barplot, pie chart, and grouped summary of site counts).

```bash
Rscript figure1/new_genes_plots.R    # -> res/plots/new_genes_stacked.pdf, new_genes_pie.pdf, new_genes_grouped.pdf
```

**Ajax/Teucer analysis** -- Metaprofile of CTCF site enrichment around the Ajax/Teucer sequence families, compared against a permuted null distribution.

```bash
Rscript figure1/ajax_teucer_analysis.R   # -> res/plots/ajax_teucer_analysis.pdf
```

#### Figure 2

**Acrocentric short arms** -- Karyoplot of the acrocentric short-arm regions (chr13/14/15/21/22), annotated with CenSat compartments and CTCF density instead of a standard cytoband track.

```bash
Rscript figure2/acrocentric_short_arms.R   # -> res/plots/acrocentric_regions.pdf
```

**Acrocentric CTCF sharing** -- Bar plot of CTCF site distribution across acrocentric short arms and a Circos plot showing CTCF site sharing between arms.

```bash
Rscript figure2/acrocentric_ctcf_sharing.R   # -> res/plots/ctcf_acro_arm_distribution.pdf, acro_circos.pdf
```

**CTCF-SDI correlation** -- Relationship between segmental duplication identity (SDI) and CTCF site density.

```bash
Rscript figure2/ctcf_sdi_correlation.R   # -> res/plots/acro_ctcf_sdi_correlation.pdf
```

**PHR CTCF coverage** -- Metaplot of CTCF occupancy across proximal heterochromatic regions (PHRs), compared against a permuted null distribution.

```bash
Rscript figure2/ctcf_phr_coverage.R   # -> res/plots/phr_ctcf_coverage.pdf
```

**rDNA plots** -- Zoomed karyoplots of CTCF and CenSat density around the chr15 acrocentric rDNA arrays.

```bash
Rscript figure2/plot_rdna.R   # -> res/plots/chr15_acro_zoom1.pdf, chr15_acro_zoom2.pdf, chr15_acro_zoom3.pdf
```

#### Figure 3

**PAF visualization** -- Alignment plots (miroPlot/SVbyEye style) of primate genomes against CHM13 target regions, built from `halLiftover`/PAF alignments.

```bash
Rscript figure3/paf_visualization.R   # -> res/plots/primate_paf_alignment.pdf
```

**CenSat conservation** -- Ridgeline plot of CTCF motif conservation scores across CenSat compartments and primates.

```bash
Rscript figure3/censat_conserv.R   # -> res/plots/primate_censat_conserv.svg
```

**CenSat conservation per chromosome** -- Same conservation analysis broken out per chromosome.

```bash
Rscript figure3/censat_conserv_per_chrom.R   # -> res/plots/primate_censat_conserv_per_chrom.svg
```

**CenSat density across primates** -- CTCF density across CenSat compartments compared across all primate assemblies.

```bash
Rscript figure3/censat_density_plots.R   # -> res/plots/density_across_primates.svg
```

**Conservation beeswarm** -- Beeswarm plot of the full distribution of CTCF motif conservation scores.

```bash
Rscript figure3/conservation_beeswarm.R   # -> res/plots/conservation_beeswarm.png
```

**Lineage-specific segmental duplications** -- Violin plot of CTCF site conservation stratified by lineage-specific segmental duplications (SDs).

```bash
Rscript figure3/lineage_specific_SDs.R   # -> res/plots/primate_SD_violin.pdf
```

**Spacer/TE analysis** -- Stacked and grouped plots of CTCF sites relative to spacer and transposable element (TE) annotations.

```bash
Rscript figure3/spacer_te_analysis.R   # -> res/plots/primate_spacer_stacked.pdf, primate_spacer_te_grouped.pdf
```

**Subterminal analysis** -- Combined CTCF occupancy and methylation metaprofiles across subterminal windows for each primate.

```bash
Rscript figure3/subterm_analysis.R   # -> res/plots/subterm_profiles.pdf
```

#### Additional analyses

**CTCF conservation metaplot** -- Heatmap of CTCF motif conservation across primates.

```bash
Rscript ctcf_conservation_metaplot.R   # -> res/plots/ctcf_conservation_heatmap.pdf
```

**DiMeLo-seq supplementary figure** -- Supplementary methylation plots (Teucer, C+T, rDNA/rDNA-promoter regions) from DiMeLo-seq data.

```bash
Rscript dimelo_supp.R   # -> res/plots/dimelo_teucer.pdf, dimelo_ct.pdf, dimelo_rdna_promoters.pdf, dimelo_rdna.pdf
```

#### Permutation tests

`src/permutation_tests/prepare_*.R` scripts prepare BED files for GAT (Genomic Association Tester) permutation testing (see [Notes](#notes) below for the GAT environment setup), and `gat_analysis.sh` runs the enrichment tests:

```bash
Rscript permutation_tests/prepare_permutation_data.R
Rscript permutation_tests/prepare_gsat_permutation_data.R
Rscript permutation_tests/prepare_phr_permutation_data.R
Rscript permutation_tests/prepare_primate_data.R
Rscript permutation_tests/prepare_teucer_permutation_data.R
bash permutation_tests/gat_analysis.sh
```

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- [HAL](https://github.com/ComparativeGenomicsToolkit/hal)
- [MEME Suite](https://meme-suite.org/meme/)

### Environment setup

```bash
conda env create -f environment.yml
conda activate ctcf
```

### Configuration

Edit `config.yml` to set the paths for your system.

## Notes

**GAT (Genomic Association Tester) setup** -- used by `src/permutation_tests/gat_analysis.sh` for permutation testing; requires its own environment due to pinned dependency versions:

```bash
conda create -n gat python=3.9 "cython<3" "numpy<1.24" -y
conda activate gat
CC=gcc pip install gat --no-build-isolation
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).