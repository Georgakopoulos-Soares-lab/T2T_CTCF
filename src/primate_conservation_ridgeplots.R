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

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results
fimo_res <- fs::path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data <- fs::path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
plots <- fs::path(config$plots) # Path to store any generated plots

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

# Store the BigWig file containg conservation scores among primates as a pointer for fast access
conserv_bw <- BigWigFile(fs::path(data, "primates", "approach2.scores_v0.3.bw"))

# Sort CHM13 FIMO results by genomic coordinates to improve the I/O speed of the 
# conservation BigWig
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
    BPPARAM = MulticoreParam(workers = 24))

# Add the conservation score to the GRanges object and convert to a dataframe
chm13_gr$mean_conservation <- unlist(conserv_scores)
chm13_df <- as.data.frame(chm13_gr)

# Get chromosome information (like chrom sizes)
chrom_info <- seqinfo(conserv_bw)

# Find the mean conservation across each chromosome in parallel
chr_conserv_means <- bplapply(
    seqnames(chrom_info), 
    function(chr) {
        gr <- import(conserv_bw, which = GRanges(chr, IRanges(1, seqlengths(chrom_info)[chr])))
        mean(gr$score, na.rm = TRUE)
    }, 
    BPPARAM = MulticoreParam(workers = 24))

# Find the genome-wide conservation mean
genome_mean <- mean(unlist(chr_conserv_means))

# Define a color palette for plotting
motif_palette <- c("MA0139.2" = "#A3ABC0", "MA1929.2" = "#424D6A", "MA1930.2" = "#0E1933")

# Create a ridgeline plot showing the distribution of CTCF site conservation scores
conserv_ridge <- chm13_df %>%
  ggplot(aes(x = mean_conservation, y = motif_id, fill = motif_id)) +
    geom_density_ridges2(
      alpha = 1,
      scale = 1.2,
      from = 0, 
      to = 1,
      color = "black",
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = genome_mean,
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
    scale_fill_manual(
      values = motif_palette
    ) +
    theme_classic() +
    labs(
      x = "Mean motif conservation",
      y = NULL,
      fill = "") +
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
  fs::path(plots, "conservation_ridgeline.pdf"),
  plot = conserv_ridge,
  width = 3,
  height = 2)