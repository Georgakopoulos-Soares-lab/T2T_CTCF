# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

gat_dir <- fs::path(res_path, "misc", "gat")

# Read RepeatMasker transposable element BED file
# Format it and then convert to GRanges
rmsk_df <- read_tsv(
  fs::path(data_path, "chm13", "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"),
  show_col_types = FALSE,
  col_names = c(
    "seqnames", "start", "end",
    "name", "score", "strand",
    "class", "subclass",
    "perc_subs", "index"
  )
)

# Filter to keep only Teucer satellites and convert to GRanges
teucer_sats <- rmsk_df |>
  dplyr::filter(grepl("teucer", name))

teucer_sats_gr <- makeGRangesFromDataFrame(
  teucer_sats,
  keep.extra.columns = FALSE
)

strand(teucer_sats_gr) <- "*"

rtracklayer::export(
  teucer_sats_gr,
  fs::path(gat_dir, "teucer.bed"),
  format = "BED"
)