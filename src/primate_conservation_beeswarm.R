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
    library("ggbeeswarm", quietly = TRUE, verbose = FALSE)
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

# Define a color palette for plotting
motif_palette <- c("MA0139.2" = "#A3ABC0", "MA1929.2" = "#424D6A", "MA1930.2" = "#0E1933")

# Create a beeswarm plot showing the full distribution of CTCF conservation scores
conserv_beeswarm <- chm13_df |>
  ggplot(
    aes(x = mean_conservation, y = motif_id, color = motif_id)
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
  geom_quasirandom(
    orientation = "y",
    size = 0.3,
    alpha = 0.4
  ) +
  geom_boxplot(
    aes(group = motif_id),
    width = 0.15,
    fill = "white",
    alpha = 0.7,
    color = "black",
    linewidth = 0.3,
    outlier.shape = NA
  ) +
  scale_color_manual(
    values = motif_palette
  ) +
  theme_classic() +
  labs(
    x = "Mean conservation",
    y = NULL,
    color = ""
  ) +
  theme(
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
    legend.position = "none",
    plot.background = NULL,
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 8),
  )

ggsave(
  fs::path(plots, "conservation_beeswarm.png"),
  plot = conserv_beeswarm,
  width = 3,
  height = 2,
  dpi = 1200,
  bg = "transparent"
)