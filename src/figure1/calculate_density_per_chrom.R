# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
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

# Load CHM13 chromosome sizes
chm13_chrom_sizes <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0.chrom.sizes"),
  show_col_types = FALSE,
  col_names = c("chr", "size")
) |>
  dplyr::filter(chr != "chrM")

# Add chromosome start (at 0) for each chromosome
chm13_chrom_sizes <- chm13_chrom_sizes |>
  dplyr::mutate(start = 1, .after = chr)

chromosome_densities <- chm13_gr |>
    as.data.frame() |>
    group_by(seqnames) |>
    count() |>
    rename(chr = seqnames) |>
    left_join(chm13_chrom_sizes, by = "chr") |>
    mutate(density = round((n / size) * 1e6)) |>
    arrange(desc(density)) |>
    rename(
        Chromosome = chr,
        `CTCF sites` = n,
        Length = size,
        `Density (per Mb)` = density
    ) |>
    dplyr::select(-c(start))

write.csv(
    chromosome_densities,
    fs::path(res_path, "misc", "chm13_chromosome_densities.csv"),
    row.names = FALSE
)
