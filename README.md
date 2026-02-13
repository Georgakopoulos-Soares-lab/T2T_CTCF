# T2T_CTCF

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
| **Diploid human** | HG002 v1.1 | T2T diploid genome (maternal + paternal haplotypes) |
| **Human pangenome** | HPRC Release 2 | 40+ samples, each with maternal and paternal assemblies |
| **Great apes** | *Pan troglodytes* | Chimpanzee (mPanTro3 v2.0) |
| | *Pan paniscus* | Bonobo (mPanPan1 v2.0) |
| | *Gorilla gorilla* | Western gorilla (mGorGor1 v2.1) |
| | *Pongo abelii* | Bornean orangutan (mPonAbe1 v2.0) |
| | *Pongo pygmaeus* | Sumatran orangutan (mPonPyg2 v2.0) |
| | *Symphalangus syndactylus* | Siamang (mSymSyn1 v2.1) |

## Repository structure

```
T2T_CTCF/
├── config.yml                # Path configuration (edit before running)
├── environment.yml           # Conda environment specification
├── links.txt                 # AWS S3 links for human pangenome downloads
├── src/
│   ├── download_chm13_data.sh      # Download CHM13 genome + annotations
│   ├── download_hg38_data.sh       # Download GRCh38.p14 genome
│   ├── download_hg002_data.sh      # Download HG002 diploid assemblies
│   ├── download_hprc_data.sh       # Download human pangenome assemblies
│   ├── download_primate_data.sh    # Download T2T great ape genomes
│   ├── download_msa_data.sh        # Download HAL multiple sequence alignments
│   ├── generate_chrom_sizes.sh     # Index FASTAs and create .chrom.sizes files
│   ├── FIMO.smk                    # Snakemake workflow for FIMO motif scanning
│   ├── FIMO.sh                     # FIMO runner script (called by Snakemake)
│   ├── create_hg38_chm13_beds.R    # Convert FIMO results to BED format
│   ├── hal_liftover.sh             # Liftover CHM13 CTCF coords to hg38 via HAL
│   ├── venn_diagram.R              # Venn diagrams comparing CHM13 vs hg38 sites
│   ├── karyoplot.R                 # Karyotype density plot of CTCF across CHM13
│   └── phastCons.smk               # Snakemake workflow for conservation scoring
├── data/                     # Input data (not tracked in git)
└── res/                      # Results (not tracked in git)
```

## Pipeline

### 1. Data download

Download genome assemblies and annotations from their respective repositories:

```bash
cd src
bash download_chm13_data.sh
bash download_hg38_data.sh
bash download_hg002_data.sh
bash download_hprc_data.sh
bash download_primate_data.sh
bash download_msa_data.sh
```

**Data sources:**
- CHM13 & HG002: [T2T Consortium](https://github.com/marbl/CHM13) (AWS S3)
- hg38: [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)
- Human pangenome: [HPRC](https://humanpangenome.org/) (AWS S3)
- Great apes: [NCBI](https://www.ncbi.nlm.nih.gov/) & [GenomeArk](https://genomeark.github.io/)
- Multiple sequence alignments: 8-way T2T ape CACTUS alignment in HAL format

### 2. Sequence preprocessing

Generate FASTA indices and chromosome size files for all downloaded genomes:

```bash
bash generate_chrom_sizes.sh
```

### 3. CTCF motif scanning (FIMO)

The Snakemake workflow auto-detects all FASTA files under `data/` and all motif files under `data/motifs/`, then:

1. Creates 1st-order Markov background models per genome (`fasta-get-markov`)
2. Scans each genome for each CTCF motif using FIMO (p-value threshold: 1e-6)

```bash
cd src
snakemake -s FIMO.smk --cores <N>
```

Results are written as TSV files to `res/fimo_res/<genome>/<motif>.tsv`.

### 4. BED conversion and coordinate liftover

Convert FIMO results for CHM13 and hg38 to BED format and lift CHM13 coordinates to hg38 using the HAL alignment:

```bash
Rscript create_hg38_chm13_beds.R
bash hal_liftover.sh
```

The liftover uses `halLiftover` to map CHM13 (hs1) CTCF site coordinates onto hg38 via the 8-way T2T ape CACTUS alignment.

### 5. Comparative analysis and visualization

**Karyotype plot** -- Visualize CTCF binding site density (per 1 Mb window) across the CHM13 karyotype, overlaid with cytobands.

```bash
Rscript karyoplot.R          # -> res/plots/karyotype_plot.pdf
```

**Venn diagrams** -- Compare CTCF sites between hg38 and T2T-CHM13 (via liftover). Sites are classified as shared (overlapping within the lifted coordinates) or assembly-specific. Regions that fail liftover are counted as CHM13-unique.

```bash
Rscript venn_diagram.R       # -> res/plots/all_liftover_venns.pdf
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

Edit `config.yml` to set the paths for your system:

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
