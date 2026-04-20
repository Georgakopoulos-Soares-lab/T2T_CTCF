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
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("BiocParallel", quietly = TRUE, verbose = FALSE)
    library("ggridges", quietly = TRUE, verbose = FALSE)
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
    select(-c(
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

  # Make the results strand-agnostic
  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Store the BigWig file containg conservation
# scores among primates as a pointer for fast access
conserv_bw <- BigWigFile(fs::path(data, "primates", "approach2.scores_v0.3.bw"))

# Sort CHM13 FIMO results by genomic coordinates
# to improve the I/O speed of the conservation BigWig
chm13_gr <- sort(chm13_gr)

# Split by chromosome for parallel processing
by_chr <- split(chm13_gr, seqnames(chm13_gr))

# Collect conservation scores for each CTCF site
# For each site, we report the mean score across the site's length
conserv_scores <- bplapply(
  by_chr,
  function(chr_gr) {
    conserv_summary <- summary(conserv_bw, which = chr_gr, type = "mean")
    vapply(conserv_summary, function(x) x$score[1], numeric(1))
  },
  BPPARAM = SnowParam(workers = 6)
)

# Add the conservation score to the GRanges object and convert to a dataframe
chm13_gr$mean_conservation <- unlist(conserv_scores)
chm13_df <- as.data.frame(chm13_gr)

# This is obtained from UCSC's bigWigInfo
genome_mean_conservation <- 0.224605

# A mapping that goes from abbreviations to full names (used for plotting)
censat_mapping <- c(
  "ct" = "Centromeric\ntransitions",
  "censat" = "Other satellites",
  "gsat" = "γSat",
  "bsat" = "βSat",
  "hsat3" = "HSat3",
  "hsat2" = "HSat2",
  "hsat1b" = "HSat1B",
  "hsat1a" = "HSat1A",
  "mon" = "αSat monomeric",
  "dhor" = "αSat divergent",
  "hor" = "αSat inactive HOR",
  "hor_l" = "αSat active HOR",
  "rdna" = "rDNA array"
)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# If an HOR has L at the end of the name then
# this indicates that this is the active alpha-satellite HOR
chm13_censat_gr$type[chm13_censat_gr$type == "hor" & grepl("L\\)$", chm13_censat_gr$name)] <- "hor_l"

# Drop unnecessary collumns
mcols(chm13_censat_gr)[, c("score", "itemRgb", "thick")] <- NULL

# Find overlaps between CTCF sites and CenSat regions
hits <- findOverlaps(chm13_gr, chm13_censat_gr, ignore.strand = TRUE)

# Initialize the censat_type column as NA for all CTCF sites
chm13_gr$censat_type <- NA

# For each CTCF site that overlaps a censat region, assign the type
overlap_types <- split(
  chm13_censat_gr$type[subjectHits(hits)],
  queryHits(hits)
)

# Assign back to the CTCF sites and
# map the short codes to full CenSat names
collapsed <- vapply(overlap_types, paste, character(1), collapse = ";")

chm13_gr$censat_type <- NA_character_
chm13_gr$censat_type[as.integer(names(collapsed))] <- collapsed

# Convert the GRanges object to a dataframe for easier plotting
ctcf_sites_df <- as.data.frame(chm13_gr)

# Add the full names to each CenSat type and drop sites that do not overlap
# with CenSat types
ctcf_sites_df$censat_type <- censat_mapping[ctcf_sites_df$censat_type]
ctcf_sites_df <- ctcf_sites_df |> filter(censat_type != "NA")

# Count how many CTCF sites exist per CenSat type
ctcf_counts <- ctcf_sites_df |>
  dplyr::mutate(
    censat_type = factor(censat_type, levels = unname(censat_mapping))
  ) |>
  dplyr::count(censat_type, .drop = FALSE) |>
  tidyr::complete(censat_type = unname(censat_mapping), fill = list(n = 0)) |>
  dplyr::mutate(label = paste0("n=", n))

# Create a ridgeline plot showing the density of points based on conservation scores
censat_conserv <- ggplot(
  ctcf_sites_df,
  aes(
    x = mean_conservation,
    y = censat_type
  )
) +
  geom_vline(
    xintercept = genome_mean_conservation,
    linetype = "dashed",
    linewidth = 0.75,
    color = "#61573D"
  ) +
  geom_vline(
    xintercept = 0.75,
    linetype = "dashed",
    linewidth = 0.75,
    color = "#e1cc97"
  ) +
  geom_density_ridges(
    fill = "#727C94",
    color = "black",
    alpha = 1,
    scale = 1,
    from = 0,
    to = 1,
    linewidth = 0.75,
    rel_min_height = 0.01
  ) +
  geom_text(
    data = ctcf_counts,
    aes(x = 0.75, y = censat_type, label = label),
    hjust = 0, vjust = -0.6,
    size = 2,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Mean conservation",
    y = NULL,
    title = NULL
  ) +
  scale_y_discrete(limits = c(
    "rDNA array", "Other satellites", "Centromeric\ntransitions",
    "HSat3", "HSat2", "HSat1B", "HSat1A", "γSat", "βSat",
    "αSat divergent", "αSat monomeric", "αSat inactive HOR",
    "αSat active HOR"
  ),
  expand = expansion(add = c(0.1, 1))
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  fs::path(plots, "primate_censat_conserv.svg"),
  censat_conserv,
  width = 2.5,
  height = 4
)