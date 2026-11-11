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

# Load the new satellite annotations
new_sats <- rtracklayer::import(
  fs::path(data_path, "chm13", "chm13v2.0_new-satellites_2022DEC.bed"),
  format = "BED"
)

# Filter to keep ajax satellites
ajax_gr <- new_sats[
  new_sats$name == "ajax"
]

# Load teucer data into a dataframe
teucer_df <- readr::read_tsv(
  fs::path(data_path, "chm13", "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"),
  show_col_types = FALSE,
  col_names = c(
   "seqnames", "start", "end",
   "name", "score", "strand",
   "class", "subclass",
   "itemRgb", "blockCount",
   "blockSizes", "blockStarts",
   "RM index", "orig_out"
  )
) |>
  dplyr::select(
    seqnames, start, end, name, strand, class, subclass
  ) |>
  dplyr::filter(
    grepl("teucer", name)
  )

# Convert teucer dataframe into a GRanges object
teucer_gr <- makeGRangesFromDataFrame(
    teucer_df,
    keep.extra.columns = TRUE
)

# Keep only teucer ranges that are next to the ajax satellites
teucer_gr <- subsetByOverlaps(teucer_gr, ajax_gr, maxgap = 100, ignore.strand = TRUE)

teucer_gr$overlaps <- countOverlaps(teucer_gr, chm13_gr, ignore.strand = TRUE)
teucer_gr$density <- (teucer_gr$overlaps / width(teucer_gr)) * 1e6

# Calculate teucer stats per chromosome
teucer_chr_density <- as.data.frame(teucer_gr) |>
  group_by(seqnames) |>
  summarise(
    n_teucer = n(),
    teucer_total_bp = sum(width),
    teucer_ctcf_overlaps = sum(overlaps),
    teucer_density = round(teucer_ctcf_overlaps / teucer_total_bp * 1e6),
    .groups = "drop"
  ) |> 
  arrange(teucer_density) |>
  filter(seqnames %in% unique(seqnames(ajax_gr)))

# Calculate ajax stats per chromosome
ajax_chr_stats <- as.data.frame(ajax_gr) |>
  group_by(seqnames) |>
  summarise(
    n_ajax = n(),
    ajax_total_bp = sum(width),
  )

# Create a table with the combined teucer and ajax sats
teucer_table <- teucer_chr_density |>
    dplyr::left_join(ajax_chr_stats, by = "seqnames") |>
    dplyr::rename(
        Chromosome = seqnames,
        `Number of teucer arrays` = n_teucer,
        `Total bp assigned to teucer arrays` = teucer_total_bp,
        `CTCF sites in teucer arrays` = teucer_ctcf_overlaps,
        `Teucer array CTCF density (per MB)` = teucer_density,
        `Number of ajax arrays` = n_ajax,
        `Total bp assigned to ajax arrays` = ajax_total_bp
    )

# Save table
write_csv(
    teucer_table,
    fs::path(res_path, "misc", "teucer_stats.csv")
)