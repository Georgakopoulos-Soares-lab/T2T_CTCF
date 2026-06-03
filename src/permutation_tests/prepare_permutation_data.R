# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
    library("plyranges", quietly = TRUE, verbose = FALSE)
    library("Biostrings", quietly = TRUE, verbose = FALSE)
    library("Rsamtools", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

# Create a dir to store the data that will be used in GAT
gat_dir <- fs::path(res_path, "misc", "gat")
fs::dir_create(gat_dir)

load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "# ", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE,
  ) |>
    dplyr::rename(name = motif_id) |>
    dplyr::select(-c(
      motif_alt_id, score,
      `p-value`, `q-value`, matched_sequence
    )) # Remove unnecessary columns
  )

  # Create one complete dataframe with the results
  res <- dplyr::bind_rows(ctcf_dfs)

  # Convert the dataframe into a GRanges
  # object for downstream interval manipulation
  fimo_gr <- GenomicRanges::makeGRangesFromDataFrame(
    res,
    seqnames.field = "sequence_name",
    start.field = "start",
    end.field = "stop",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$name <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define active alpha satellites based on CenSat color
chm13_censat_gr$name[
  chm13_censat_gr$name == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Finally, keep only the type as a metadata column
chm13_censat_gr <- chm13_censat_gr[, "name"]

# Load the chromosome sizes for the CHM13 genome
chm13_chrom_sizes <- read_tsv(
  fs::path(
    data, "chm13", "chm13v2.0.chrom.sizes"
  ),
  show_col_types = FALSE,
  col_names = c("chrom", "size")
) |>
  dplyr::mutate(start = 0, workspace = "ws") |>
  dplyr::filter(chrom != "chrM") |>
  dplyr::select(chrom, start, size, workspace)

# Create a named vector that will be used to tile the genome
chm13_lengths <- setNames(
  chm13_chrom_sizes$size,
  chm13_chrom_sizes$chrom
)

# Tile the genome in 100 Kb windows
windows <- GenomicRanges::tileGenome(
  chm13_lengths,
  tilewidth = 100e3,
  cut.last.tile.in.chrom = TRUE
)

# Define the path where the T2T-CHM13v2.0 genome FASTA is store
fasta_file <- fs::path(data, "chm13", "chm13v2.0.fa.gz")

# Provide the FASTA file as a pointer for R to use
# Make sure the genome file is indexed using Samtools
fa <- Rsamtools::FaFile(fasta_file)

# Get the underlying sequence for each window
window_sequences <- Biostrings::getSeq(fa, windows)

# Calculate GC content (counts / region's length)
# and add it as metadata to the windows GRanges
gc_count <- as.integer(
  Biostrings::letterFrequency(
    window_sequences, "GC"
  )
)

windows <- windows |>
  mutate(
    gc_count = gc_count,
    gc_content = gc_count / width(windows)
  )

# Define isochores to be used with GAT
windows <- windows |>
  mutate(
    name = cut(
      gc_content,
      breaks = quantile(gc_content, probs = seq(0, 1, 1 / 6), na.rm = TRUE),
      labels = c("L1", "L2", "M1", "M2", "H1", "H2"),
      include.lowest = TRUE
    )
  ) |>
  select(-c(gc_content, gc_count))

# Export the files needed for the GAT
# Isochores are regions with similar GC%
rtracklayer::export(
  windows,
  fs::path(gat_dir, "isochores.bed"),
  format = "BED"
)

# CTCF are used as segments
rtracklayer::export(
  chm13_gr,
  fs::path(gat_dir, "ctcf.bed"),
  format = "BED"
)

# CenSat is used as annotations
rtracklayer::export(
  chm13_censat_gr,
  fs::path(gat_dir, "censat.bed"),
  format = "BED"
)

# Workspace is defined as the whole genome, where we can put CTCF sites
chm13_chrom_sizes |>
  write_tsv(
    fs::path(gat_dir, "workspace.bed"),
    col_names = FALSE
  )