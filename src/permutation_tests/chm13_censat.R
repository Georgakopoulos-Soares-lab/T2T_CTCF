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
    library("regioneR", quietly = TRUE, verbose = FALSE)
    library("plyranges", quietly = TRUE, verbose = FALSE)
    library("Biostrings", quietly = TRUE, verbose = FALSE)
    library("Rsamtools", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "# ", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE,
  ) |>
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

  GenomicRanges::strand(fimo_gr) <- "*"

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
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Finally, keep only the type as a metadata column
chm13_censat_gr <- chm13_censat_gr[, "type"]

# Load the chromosome sizes for the CHM13 genome
chm13_chrom_sizes <- read_tsv(
  fs::path(
    data, "chm13", "chm13v2.0.chrom.sizes"
  ),
  show_col_types = FALSE,
  col_names = c("chrom", "size")
)

# Create a named vector that will be used to tile the genome
chm13_lengths <- setNames(
  chm13_chrom_sizes$size,
  chm13_chrom_sizes$chrom
)

# Tile the genome in 1 Kb windows
windows <- GenomicRanges::tileGenome(
  chm13_lengths,
  tilewidth = 1000,
  cut.last.tile.in.chrom = TRUE
)

# Define the path where the T2T-CHM13v2.0 genome FASTA is store
fasta_file <- fs::path(data, "chm13", "chm13v2.0.fa.gz")

# Provide the FASTA file as a pointer for R to use
# Make sure the genome file is indexed using Samtools
fa <- Rsamtools::FaFile(fasta_file)

# Get the underlying sequence for each window
window_sequences <- Biostrings::getSeq(fa, windows)

# Calculate CpG and GC content (counts / region's length)
# and add it as metadata to the windows GRanges
cpg_counts <- Biostrings::vcountPattern(
  "CG",
  window_sequences
)

gc_counts <- as.integer(
  Biostrings::letterFrequency(
    window_sequences, "GC"
  )
)

windows <- windows |>
  mutate(
    cpg_count = cpg_counts,
    cpg_content = cpg_counts / width(windows),
    gc_count = gc_counts,
    gc_content = gc_counts / width(windows),
    idx = seq_along(windows)
  )

# For each window find, which CenSat type it overlaps the most and assign
# a label
labels <- windows |>
  join_overlap_intersect(chm13_censat_gr) |>
  as_tibble() |>
  group_by(idx) |>
  slice_max(width, n = 1, with_ties = FALSE) |>
  ungroup() |>
  dplyr::select(idx, type)

# Add this CenSat label back to the windows GRanges object
windows$type <- "non_satellite"
windows$type[labels$idx] <- as.character(labels$type)

# Create GC content bins that increment by 1%
gc_bins  <- quantile(
  mcols(windows)$gc_content,
  probs = seq(0, 1, 0.01)
)

# Create CpG content bins that increment by 10%
cpg_bins <- quantile(
  mcols(windows)$cpg_content,
  probs = seq(0, 1, 0.1)
)

# Create composite bins which contain both the CpG and GC bin information
mcols(windows)$bin <- paste(
  findInterval(
    mcols(windows)$gc_content, gc_bins, all.inside = TRUE
  ),
  findInterval(
    mcols(windows)$cpg_content, cpg_bins, all.inside = TRUE
  ),
  sep = "_"
)

# Create a list with the available CenSat types
censat_types <- unique(chm13_censat_gr$type)

# Create a list with the available CTCF motif IDs
ctcf_motifs <- unique(chm13_gr$motif_id)

RNGkind("L'Ecuyer-CMRG") # Reproducible streams across parallel workers
set.seed(42) # Set seed so that each permutation test result is reproducible
results_list <- list() # Define a list where we will store the results

# Iterate over all CenSat types
for (censat_type in censat_types) {

  # Filter to only the windows that contain the
  # CenSat category we analyze
  region_windows <- windows |>
    filter(type == censat_type)

  # Find which bins are needed for the permutation and keep them as a pool
  acceptable_bins <- unique(region_windows$bin)

  # Create a pool with the acceptable regions where the random shuffling can occur
  pool <- windows |>
    filter(bin %in% acceptable_bins)

  # Filter only the CenSat GRanges to only the regions that are analyzed
  region_gr <- chm13_censat_gr[chm13_censat_gr$type == censat_type]

  # Iterate over CTCF motifs
  for (motif in ctcf_motifs) {
    message(sprintf("Running: %s x %s", censat_type, motif))
    motif_gr <- chm13_gr[chm13_gr$motif_id == motif]
    suppressMessages(suppressWarnings({
      pt <- permTest(
        A = motif_gr, # Motif FIMO res
        B = region_gr, # Regions assigned to one CenSat type
        ntimes = 10000, # Number of permutations
        randomize.function = randomizeRegions,
        evaluate.function = numOverlaps,
        alternative = "greater", # We only test for enrichment (not depletion)
        genome = pool,
        mask = NULL,
        allow.overlaps = TRUE,
        per.chromosome = TRUE,
        verbose = FALSE,
        # Seed is removed, so that each core
        # uses the same seed defined in the start
        mc.set.seed = FALSE,
        force.parallel = TRUE,
        mc.cores = 130
      )

      # Store the results in a dataframe
      key <- paste(censat_type, motif, sep = "__")
      results_list[[key]] <- tibble(
        censat = censat_type,
        motif  = motif,
        pval   = pt$numOverlaps$pval,
        zscore = pt$numOverlaps$zscore
      )
    }))
  }
}

# Create the final dataframe
results_df <- dplyr::bind_rows(results_list, .id = "key")

# Save results to a CSV file
# We will load these results in the code
# we use for the plotting of the densities across CenSat types
write.csv(
  results_df,
  fs::path(res_path, "misc", "chm13_permtest_results.csv"),
  row.names = FALSE
)