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
        library("EnrichedHeatmap", quietly = TRUE, verbose = FALSE)
        library("magick", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

# Path containing FIMO results for CHM13,
# HG002, T2T primates, and human pangenome asseblies
fimo_res <- fs::path(config$fimo_res)

# Path containing metadata for CHM13,
# HG002, T2T primates, and human pangenome asseblies
data <- fs::path(config$data)

# Path to store any generated plots
plots <- fs::path(config$plots)

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

# Load CHM13 results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Sort CHM13 FIMO results by genomic coordinates to improve the I/O speed of the 
# conservation BigWig
chm13_gr <- sort(chm13_gr)

# Split by chromosome for parallel processing
by_chr <- split(chm13_gr, seqnames(chm13_gr))

# Define a function that adds a conservation score to each CTCF site
# The score of each site is the mean conservation score of the region
annotate_sites <- function(
    chm13_gr,
    split_gr,
    annotation_bw_path
) {
    # Ensure the BigWig path is a plain character string
    annotation_bw_path <- as.character(annotation_bw_path)

    # Parallel-iterate over each chromosome's GRanges in split_gr.
    # Each worker independently opens the BigWig file and computes
    # the mean signal score for every genomic range in its chunk.
    scores <- bplapply(
        split_gr,
        function(chr_gr, bw_path) {
            # Open a handle to the BigWig file within the worker process
            annotation_bw <- BigWigFile(bw_path)

            # Compute the mean signal across each range in chr_gr.
            # Returns a list of single-row DataFrames, one per range.
            summary <- rtracklayer::summary(annotation_bw, which = chr_gr, type = "mean")

            # Extract the numeric score from each summary DataFrame,
            # collapsing the list into a single numeric vector
            vapply(summary, function(x) x$score[1], numeric(1))
        },
        bw_path = annotation_bw_path,
        BPPARAM = SnowParam(workers = 24)
    )

    scores
}

# Add conservation scores to the GRanges object
chm13_gr$mean_conservation <- unlist(annotate_sites(
    chm13_gr,
    by_chr,
    fs::path(data, "primates", "approach2.scores_v0.3.bw")
))

# Extract conserved CTCF sites (defined as those with a score greater or equal to 75%)
conserved_ctcf_sites <- chm13_gr[mcols(chm13_gr)$mean_conservation >= 0.75]

# Extend motifs by 50kb
extended_conserved_ctcf_sites <- conserved_ctcf_sites + 50000

# Load conservation scores for the extented CTCF regions
conserv_scores <- import(fs::path(data, "primates", "approach2.scores_v0.3.bw"), which = extended_conserved_ctcf_sites)

# Calculate matrix up to 10kb to see if the CTCF conservation is localized or not
conserv_mat = normalizeToMatrix(
    conserv_scores,
    conserved_ctcf_sites,
    value_column = "score", 
    extend = 10000,
    mean_mode = "w0",
    w = 50
)

# Save the metaplot to a PDF file
pdf(fs::path(plots, "ctcf_conservation_heatmap.pdf"), width = 3, height = 4)

EnrichedHeatmap(
    conserv_mat,
    name = "CTCF Conservation",
    col = c("white", "#0E1933"),
    axis_name = c("-10kb", "CTCF site", "", "10kb"),
    use_raster = TRUE,
    top_annotation = HeatmapAnnotation(
        enriched = anno_enriched(gp = gpar(col = "#0E1933", lwd = 2))
    )
)

dev.off()
