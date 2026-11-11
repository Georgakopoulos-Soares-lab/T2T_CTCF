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
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)

load_fimo_results <- function(genome, fimo_res_dir) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(fs::path(fimo_res_dir, genome), glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
  ) |>
    select(
      -c(motif_alt_id, score, `p-value`, `q-value`, matched_sequence)
    ) # Remove unnecessary columns
  )

  # Create one complete dataframe with the results
  res <- dplyr::bind_rows(ctcf_dfs)

  # Convert the dataframe into a GRanges object
  # for downstream interval manipulation
  ctcf_gr <- GenomicRanges::makeGRangesFromDataFrame(
    res,
    seqnames.field = "sequence_name",
    start.field = "start",
    end.field = "stop",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

  ctcf_gr
}

chm13_gr <- load_fimo_results("chm13", fimo_res)

# Load CHM13 CenSat annotations
chm13_censat_gr <- rtracklayer::import(
  fs::path(data_path, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define alpha satellite HOR as the region that has the associated color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Find CTCF sites overlapping CT regions
ct_gr <- chm13_censat_gr[chm13_censat_gr$type == "ct"]
ct_ctcf <- subsetByOverlaps(chm13_gr, ct_gr, ignore.strand = TRUE)

# Load CHM13 segmental duplication annotations
seg_dups <- import(
  fs::path(data_path, "chm13", "chm13v2.0_SD.bed"),
  format = "BED"
)

ct_ctcf_in_seg_dups <- subsetByOverlaps(ct_ctcf, seg_dups, ignore.strand = TRUE)

perc_ct_ctcf_in_seq_dups <- round(
    nrow(as.data.frame(ct_ctcf_in_seg_dups)) / nrow(as.data.frame(ct_ctcf)) * 100,
    2
)