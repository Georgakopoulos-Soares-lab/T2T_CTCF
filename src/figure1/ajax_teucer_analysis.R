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
    library("EnrichedHeatmap", quietly = TRUE, verbose = FALSE)
    library("tibble", quietly = TRUE, verbose = FALSE)
    library("regioneR", quietly = TRUE, verbose = FALSE)
    library("BiocParallel", quietly = TRUE, verbose = FALSE)
    library("matrixStats", quietly = TRUE, verbose = FALSE)
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
ajax_sats <- new_sats[
  new_sats$name == "ajax"
]

# Create a matrix using EnrichedHeatmap to
# visualize CTCF density within ajax satellites
ajax_mat <- EnrichedHeatmap::normalizeToMatrix(
  signal = chm13_gr,
  target = ajax_sats,
  extend = 1500,
  w = 15,
  mean_mode = "absolute",
  value_column = NULL,
  background = 0,
  target_ratio = 0.4,
  verbose = FALSE,
  smooth = FALSE
)

# Load GC isochores 10 Kb non-overlapping tiles
isochores <- rtracklayer::import(
  fs::path(res_path, "misc", "gat", "isochores.bed")
)

N_ITER <- 10 # number of permutations
SEED <- 42 # random seed so that results are reprodusible
WORKERS <- 4 # number of cores to use
BPP <- MulticoreParam(workers = WORKERS)

# Create a list with the arguments that are used to create the coverage matrix
common_args <- list(
  target = ajax_sats,
  extend = 1500,
  w = 15,
  mean_mode = "absolute",
  value_column = NULL,
  background = 0,
  target_ratio = 0.4,
  verbose = FALSE,
  smooth = FALSE
)

# Helper function that assign isochores to CTCF sites
assign_isochores <- function(
  ctcf_gr, isochores
) {

  # Find with which isochore each CTCF sites overlaps with
  hits <- GenomicRanges::findOverlaps(
    ctcf_gr,
    isochores,
    select = "first"
  )

  # Assign the isochore label back to the original
  # GRanges object
  ctcf_gr$name <- isochores$name[hits]

  # Return the updated object
  ctcf_gr
}

# Assign isochores to short arm CTCF sites present in
# one or multiple arms
chm13_gr_with_isochores <- assign_isochores(
  chm13_gr,
  isochores
)

# Helper function that generates a GRangesList of GC matched shuffled CTCF sites
shuffle_within_isochores <- function(
  regions,
  isochores,
  n_iter,
  seed,
  BPPARAM = BiocParallel::bpparam()
) {

  # Extract the availabe isochore labels
  isochore_labels <- unique(
    GenomicRanges::mcols(regions)$name
  )

  # Set a seed that will be reproducible across cores
  BiocParallel::bpRNGseed(BPPARAM) <- seed

  # Generate shuffled sets in parallel
  shuffles <- BiocParallel::bplapply(
    seq_len(n_iter), # go over iterations
    function(i) {
      # For each isochore label
      per_label <- lapply(isochore_labels, function(label) {
        # extract CTCF sites that match the label
        ctcf_subset <- regions[
          GenomicRanges::mcols(regions)$name == label
        ]

        # extract isochores that match the label
        isochore_subset <- isochores[
          GenomicRanges::mcols(isochores)$name == label
        ]

        # use regioneR to generate null sets matching isochores
        # We allow overlapping regions and the permutation is
        # performed with no chromosome constriction
        regioneR::randomizeRegions(
          ctcf_subset,
          genome = isochore_subset,
          mask = NA,
          allow.overlaps = TRUE,
          per.chromosome = FALSE
        )
      })

      # Create the shuffled set
      shuffled_gr <- do.call(c, per_label)

      # Add the iteration number to keep track
      GenomicRanges::mcols(shuffled_gr)$iter <- i

      # Return the objevt
      shuffled_gr
    },
    BPPARAM = BPPARAM
  )

  # Create the GRanges list
  GenomicRanges::GRangesList(shuffles)
}

# Generate null sets for both single arm and
# multiple arm sites
null_set <- shuffle_within_isochores(
  chm13_gr_with_isochores, isochores,
  n_iter = N_ITER, seed = SEED, BPPARAM = BPP
)

# Helper function that runs EnrichedHeatmap's normalizeToMatrix
# for a null region set and returns the mean value per column (bin)
col_means_for_signal <- function(signal_gr, common_args) {

  # Call function
  m <- do.call(
    EnrichedHeatmap::normalizeToMatrix,
    c(list(signal = signal_gr), common_args)
  )

  # Get a per-bin mean signal across all target regions
  colMeans(m)
}

# Helper function that runs the matrix generation in parallel
null_matrix <- function(null_set, common_args, BPPARAM = bpparam()) {
  rows <- BiocParallel::bplapply(
    seq_along(null_set),
    function(i) col_means_for_signal(null_set[[i]], common_args),
    BPPARAM = BPPARAM
  )
  do.call(rbind, rows)
}

# Create the associated matrices for the null sets
null_mat <- null_matrix(
  null_set,
  common_args,
  BPP
)

# Helper function that summarizes the results into a single dataframe for
# plotting. We also keep track of the null distribution's centre (mean/median)
# and a pointwise 95% interval — the 2.5th and 97.5th percentiles of the null
# at each bin — which become the grey ribbon the observed mean is judged against
summarise_with_null <- function(observed_mat, null_mat) {
  observed <- colMeans(observed_mat, na.rm = TRUE)

  tibble::tibble(
    bin = seq_len(ncol(observed_mat)),
    mean = observed,
    null_mean = colMeans(null_mat),
    null_median = matrixStats::colMedians(null_mat),
    null_q025 = matrixStats::colQuantiles(null_mat, probs = 0.025),
    null_q975 = matrixStats::colQuantiles(null_mat, probs = 0.975)
  )
}

# Create the dataframe we will use for plotting
plot_df <- dplyr::bind_rows(
  summarise_with_null(
    ajax_mat,
    null_mat
  )
)

# Extract the indices for where the normalized region starts and ends
# and for where the flanking regions end
upstream_n <- length(attr(ajax_mat, "upstream_index"))
target_idx <- attr(ajax_mat, "target_index")
target_n <- length(target_idx)
target_start <- min(target_idx)
target_end   <- max(target_idx)
n_bins <- ncol(ajax_mat)

# Create the metaprofile plot and save to a PDF
ajax_metaplot <- ggplot(
  plot_df, aes(x = bin)
) +
  geom_line(
    data = subset(plot_df, bin <= target_start),
    aes(y = mean),
    color = "#274a96", linewidth = 0.5
  ) +
  geom_line(
    data = subset(plot_df, bin >= target_start & bin <= target_end),
    aes(y = mean),
    color = "#967527", linewidth = 0.5
  ) +
  geom_line(
    data = subset(plot_df, bin >= target_end),
    aes(y = mean),
    color = "#274a96", linewidth = 0.5
  ) +
  geom_ribbon(
    aes(ymin = null_q025, ymax = null_q975),
    fill = "grey70",
    colour = NA,
    alpha = 0.35
  ) +
  geom_line(
    aes(y = null_median),
    colour = "grey50",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = c(1, upstream_n, upstream_n + target_n, n_bins),
    labels = c("-1.5 Kb", "Start", "End", "+1.5 Kb"),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(0.25),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(x = "Ajax-Teucer arrays", y = "CTCF coverage") +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    plot.background = NULL,
    legend.position = "none"
  )

ggsave(
  fs::path(res_path, "plots", "ajax_teucer_analysis.pdf"),
  ajax_metaplot,
  width = 2.5,
  height = 1.8,
  bg = "transparent"
)