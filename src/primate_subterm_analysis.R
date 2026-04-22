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
    library("magick", quietly = TRUE, verbose = FALSE)
    library("patchwork", quietly = TRUE, verbose = FALSE)
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

  # Make the object strand-agnostic
  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

# Create a function that plots CTCF coverage in subterminal heterochromatin
# regions and their spacers (just showing 30Kb flanks)
plot_subterm_profile <- function(
  species_dir,
  censat_file,
  meth_file,
  title,
  regex,
  show_x_labels = TRUE,
  show_y_labels = TRUE
) {
  # Load CTCF sites
  gr <- load_fimo_results(
    fs::path(fimo_res, "primates", species_dir)
  )

  # Load CenSat annotations
  censat_gr <- rtracklayer::import(
    fs::path(data, "primates", species_dir, censat_file),
    format = "BED"
  )

  # Load methylation data
  meth_bw <- rtracklayer::import(
    fs::path(data, "primates", species_dir, meth_file),
    format = "BigWig"
  )

  # Add a column specifying the type of each region
  censat_gr$type <- tolower(sub(regex, "", censat_gr$name))

  # Remove unused columns and keep only subterminal heterochromatin regions
  GenomicRanges::mcols(censat_gr)[, c("score", "itemRgb", "thick")] <- NULL
  subterm_gr <- censat_gr[censat_gr$type == "subterm"]

  # Create a matrix showing CTCF occupancy in 100bp windows across
  # subterminal heterochromatin regions and 30Kb flanks
  fimo_mat <- EnrichedHeatmap::normalizeToMatrix(
    gr, subterm_gr,
    extend = 30000, mean_mode = "absolute",
    w = 100, background = 0, smooth = FALSE
  )

  # Count the methylation levels in 100bp windows across
  # subterminal heterochromatin regions and 30Kb flanks
  meth_mat <- EnrichedHeatmap::normalizeToMatrix(
    meth_bw, subterm_gr,
    value_column = "score", mean_mode = "absolute",
    extend = 30000, w = 100, background = 0
  )

  # Convert both matrices into dataframes
  fimo_cov <- data.frame(
    pos = seq_len(ncol(fimo_mat)),
    signal = colMeans(fimo_mat)
  )
  meth_cov <- data.frame(
    pos = seq_len(ncol(meth_mat)),
    signal = colMeans(meth_mat)
  )

  # Find the indices describing the position of the
  # subterminal heterochromatin region
  target_idx <- attr(fimo_mat, "target_index")
  target_start <- min(target_idx) # Find the start
  target_end   <- max(target_idx) # Find the end
  x_breaks <- c(1, min(target_idx), max(target_idx), ncol(fimo_mat))

  # Line plot showing CTCF occupancy in the windows
  fimo_plot <- ggplot2::ggplot(
    fimo_cov, ggplot2::aes(pos, signal)
  ) +
    ggplot2::geom_line(
      data = subset(fimo_cov, pos <= target_start),
      color = "#9D5237", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = subset(fimo_cov, pos >= target_start & pos <= target_end),
      color = "#22233F", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = subset(fimo_cov, pos >= target_end),
      color = "#9D5237", linewidth = 0.3
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (show_y_labels) "Mean\ncoverage" else NULL,
      title = title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = if (show_y_labels) {
        ggplot2::element_text(size = 6)
      } else {
        ggplot2::element_blank()
      },
      axis.title.y = if (show_y_labels) {
        ggplot2::element_text(size = 8)
      } else {
        ggplot2::element_blank()
      },
      panel.grid.major.y = ggplot2::element_line(
        color = "gray70",
        linewidth = 0.1
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "gray70",
        linewidth = 0.1
      ),
      plot.title = ggplot2::element_text(face = "italic", size = 8)
    )

  # Line plot showing methylation in the windows
  meth_plot <- ggplot2::ggplot(
    meth_cov, ggplot2::aes(pos, signal)
  ) +
    ggplot2::geom_hline(
      yintercept = 0.3,
      linetype = "dashed",
      color = "#247cd3",
      linewidth = 0.2
    ) +
    ggplot2::geom_hline(
      yintercept = 0.7,
      linetype = "dashed",
      color = "#ca2e40",
      linewidth = 0.2
    ) +
    ggplot2::geom_line(
      data = subset(meth_cov, pos <= target_start),
      color = "#9D5237", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = subset(meth_cov, pos >= target_start & pos <= target_end),
      color = "#22233F", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      data = subset(meth_cov, pos >= target_end),
      color = "#9D5237", linewidth = 0.3
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = if (show_x_labels) {
        c("-30Kb", "start", "end", "30Kb")
      } else {
        rep("", length(x_breaks))
      }
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0.0", "0.5", "1.0")
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (show_y_labels) "Mean\nmethylation" else NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        color = "gray70",
        linewidth = 0.1
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "gray70",
        linewidth = 0.1
      ),
      axis.text.y = if (show_y_labels) {
        ggplot2::element_text(size = 6)
      } else {
        ggplot2::element_blank()
      },
      axis.title.y = if (show_y_labels) {
        ggplot2::element_text(size = 8)
      } else {
        ggplot2::element_blank()
      },
      axis.text.x = if (show_x_labels) {
        ggplot2::element_text(size = 8)
      } else {
        ggplot2::element_blank()
      }
    )

  # Combine the 2 plots into one
  combined <- (fimo_plot / meth_plot) +
    patchwork::plot_layout(heights = c(2, 1))

  list(
    plot = combined,
    plots = list(coverage = fimo_plot, methylation = meth_plot),
    data = list(coverage = fimo_cov, methylation = meth_cov),
    matrices = list(fimo = fimo_mat, meth = meth_mat),
    subterm_gr = subterm_gr
  )
}

# Create the plots for each primate
suppressWarnings({
  gorgor_plot <- plot_subterm_profile(
    "GorGor",
    "mGorGor1_v2.0_CenSat_v2.0.bed",
    "mGorGor1_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw",
    "Gorilla gorilla gorilla",
    "\\(.*",
    TRUE,
    TRUE
  )

  panpan_plot <- plot_subterm_profile(
    "PanPan",
    "mPanPan1_v2.0_CenSat_v2.0.bed",
    "mPanPan1_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw",
    "Pan paniscus",
    "\\(.*",
    FALSE,
    TRUE
  )

  pantro_plot <- plot_subterm_profile(
    "PanTro",
    "mPanTro3_v2.0_CenSat_v2.0.bed",
    "mPanTro3_v2.0.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw",
    "Pan troglodytes",
    "\\(.*",
    FALSE,
    FALSE
  )

  symsyn_plot <- plot_subterm_profile(
    "SymSyn",
    "mSymSyn1_v2.1_CenSat_v2.0.bed",
    "mSymSyn1_v2.1.ONT_Guppy6.3.8_5mC.winnowmap_modkit_5mC_CpG_v1.0.bw",
    "Symphalangus syndactylus",
    "_.*",
    TRUE,
    FALSE
  )
})

# Combine individual plots into one and save
combined_profiles <- (
  panpan_plot$plot | pantro_plot$plot
) / (
  gorgor_plot$plot | symsyn_plot$plot
)

ggsave(
  fs::path(plots, "subterm_profiles.pdf"),
  combined_profiles,
  width = 5.5,
  height = 3.5,
  bg = "transparent"
)