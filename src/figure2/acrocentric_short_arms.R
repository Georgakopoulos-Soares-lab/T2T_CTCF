# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("karyoploteR", quietly = TRUE, verbose = FALSE)
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
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

  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))
chm13_censat_gr$annotation <- tolower(stringr::str_extract(
  chm13_censat_gr$name, "(?<=\\().*?(?=\\))"
))

chm13_censat_gr$annotation[is.na(chm13_censat_gr$annotation)] <- "No annotation"

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Define a function where a pointer to BigWig files is created
point_to_bw <- function(name) {
  bw <- rtracklayer::BigWigFile(
    fs::path(
      data,
      "chm13",
      name
    )
  )

  bw
}

# Keep the long-read methylation file as a pointer
methylation_bw <- point_to_bw("chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw")

# Keep the Fiber-Seq accessibility file as a pointer
accessibility_bw <- point_to_bw("all.percent.accessible.bw")

# Long-read DiMeLo-seq histone modifications
h3k27ac_bw <- point_to_bw("h3k27ac_dimelo.bw")
h3k27me3_bw <- point_to_bw("h3k27me3_dimelo.bw")
h3k4me3_bw <- point_to_bw("h3k4me3_dimelo.bw")

# Extract the acrocentric short arm regions for plotting
acrocentric_regions <- subset(
  chm13_censat_gr,
  seqnames %in% c("chr13", "chr14", "chr15", "chr21", "chr22") & type == "hor_l"
)

# Set all starts to 1
start(acrocentric_regions) <- 1

# Create one unified region for each chromosome
acrocentric_regions <- unlist(
  range(
    split(
      acrocentric_regions,
      seqnames(acrocentric_regions)
    )
  )
)

# Extract the CenSat annotations assigned only to the
# acrocentric short arms
acrocentric_censat <- subsetByOverlaps(
  chm13_censat_gr,
  acrocentric_regions,
  ignore.strand = TRUE
)

pdf(
  fs::path(plots, "acrocentric_regions.pdf"),
  width = 8.5,
  height = 6.5
)

suppressWarnings({
  # Modify plot parameters
  pp <- karyoploteR::getDefaultPlotParams(plot.type = 1)
  pp$ideogramheight <- 15
  pp$data1inmargin <- 10
  pp$rightmargin <- 0.03
  pp$leftmargin <- 0.18
  pp$topmargin <- 0
  pp$bottommargin <- 32

  # Initialize plot with no ideogram since we will populate that
  # with CenSat annotations
  kp <- karyoploteR::plotKaryotype(
    genome = acrocentric_regions,
    plot.params = pp,
    ideogram.plotter = NULL,
    cex = 0.9,
    chromosomes = "all"
  )

  # Instead of a cytoband track plot the CenSat annotations
  kp <- karyoploteR::kpPlotRegions(
    kp,
    acrocentric_censat,
    r0 = 0,
    r1 = 1,
    col = acrocentric_censat$itemRgb,
    border = NA,
    data.panel = "ideogram",
    clipping = TRUE
  )

  # Add base numbering every 5 Mb with minor ticks every 1 Mb
  kp <- karyoploteR::kpAddBaseNumbers(
    kp,
    tick.dist = 5e6,
    minor.tick.dist = 1e6,
    tick.len = 5,
    add.units = TRUE,
    digits = 2,
    cex = 0.4
  )

  # Pre-compute CTCF density (we still need kpPlotDensity to get the binning,
  # but we plot it invisibly first, then re-plot per chromosome with proper axes)

  window <- 1e4  # window size for aggregation

  # Dummy density pass to get per-window CTCF counts
  kp <- karyoploteR::kpPlotDensity(
    kp,
    data = chm13_gr,
    col = NA,
    border = NA,
    window.size = window,
    r1 = 0.01
  )

  # Extract the density values and get the per chromosome
  # max densities
  computed <- kp$latest.plot$computed.values

  per_chrom_ctcf_max <- tapply(
    computed$density,
    as.character(GenomicRanges::seqnames(computed$windows)),
    max
  )

  # Define track heights
  ctcf_r0 <- 0.00
  ctcf_r1 <- 0.14

  acc_r0 <- 0.17
  acc_r1 <- 0.31

  meth_r0 <- 0.34
  meth_r1 <- 0.48

  h3k4me3_r0 <- 0.51
  h3k4me3_r1 <- 0.65

  h3k27ac_r0  <- 0.68
  h3k27ac_r1  <- 0.82

  h3k27me3_r0  <- 0.85
  h3k27me3_r1  <- 0.99

  # Bin the acrocentric short arm regions in length equal to the
  # specified window
  short_arm_bins <- unlist(
    GenomicRanges::tile(
      acrocentric_regions,
      width = window
    )
  )

  # Define a helper value to extract the max value of a
  # specific track per each chromosome
  bin_bigwig_scores <- function(
    bigwig, bins
  ) {
    bw_scores <- unlist(
      rtracklayer::summary(
        object = bigwig,
        which  = bins,
        size   = 1,
        type   = "mean"
      )
    )

    per_chrom_max <- tapply(
      bw_scores$score,
      as.character(GenomicRanges::seqnames(bw_scores)),
      max, na.rm = TRUE
    )

    list(bw_scores, per_chrom_max)
  }

  # Extract per chromosome Fiber-seq accessibility binned scores
  acc_list <- bin_bigwig_scores(
    accessibility_bw,
    short_arm_bins
  )

  # Extract per chromosome Nanopore methylation binned scores
  meth_list <- bin_bigwig_scores(
    methylation_bw,
    short_arm_bins
  )

  # Extract per chromosome histone modification tracks
  h3k27ac_list <- bin_bigwig_scores(
    h3k27ac_bw,
    short_arm_bins
  )

  h3k27me3_list <- bin_bigwig_scores(
    h3k27me3_bw,
    short_arm_bins
  )

  h3k4me3_list <- bin_bigwig_scores(
    h3k4me3_bw,
    short_arm_bins
  )

  # Create a helper function to plot the binned BigWig files
  plot_bigwig <- function(
    kp, scores_list, chrom,
    track_r0, track_r1, label,
    color
  ) {
    chrom_max <- scores_list[[2]][[chrom]]
    chrom_values <- scores_list[[1]][
      as.character(GenomicRanges::seqnames(scores_list[[1]])) == chrom
    ]

    chrom_values$score[is.na(chrom_values$score)] <- 0

    kp <- karyoploteR::kpAxis(
      kp,
      chromosomes = chrom,
      ymin = 0,
      ymax = ceiling(chrom_max),
      r0 = track_r0,
      r1 = track_r1,
      data.panel = 1,
      cex = 0.45,
      tick.pos = ceiling(chrom_max),
      side = 1
    )

    kp <- karyoploteR::kpArea(
      kp,
      data = chrom_values,
      y    = chrom_values$score,
      ymin = 0,
      ymax = chrom_max,
      col    = color,
      border = color,
      lwd    = 0.5,
      r0 = track_r0,
      r1 = track_r1,
      data.panel = 1,
      clipping   = TRUE
    )

    kp <- karyoploteR::kpAddLabels(
      kp,
      labels = label,
      r0 = track_r0,
      r1 = track_r1,
      data.panel = 1,
      side = "left",
      cex = 0.6,
      label.margin = 0.03
    )
  }

  # Plot the tracks for each chromosome separetely
  for (chrom in names(per_chrom_ctcf_max)) {
    chrom_ctcf_max <- per_chrom_ctcf_max[[chrom]]
    chrom_ctcf <- chm13_gr[
      GenomicRanges::seqnames(chm13_gr) == chrom
    ]

    kp <- karyoploteR::kpAxis(
      kp,
      chromosomes = chrom,
      ymin = 0,
      ymax = chrom_ctcf_max,
      r0 = ctcf_r0,
      r1 = ctcf_r1,
      data.panel = 1,
      cex = 0.45,
      tick.pos = chrom_ctcf_max,
      side = 1
    )

    kp <- karyoploteR::kpPlotDensity(
      kp,
      data = chrom_ctcf,
      window.size = window,
      ymax = chrom_ctcf_max,
      col = "#4D4D4D",
      border = "gray30",
      lwd = 0.5,
      r0 = ctcf_r0,
      r1 = ctcf_r1,
      data.panel = 1,
      clipping = TRUE
    )

    kp <- karyoploteR::kpAddLabels(
      kp,
      labels = "CTCF density",
      r0 = ctcf_r0,
      r1 = ctcf_r1,
      data.panel = 1,
      side = "left",
      cex = 0.6,
      label.margin = 0.03
    )

    kp <- plot_bigwig(
      kp, acc_list, chrom,
      acc_r0, acc_r1, "Accessible chromatins (%)",
      "#608266"
    )

    kp <- plot_bigwig(
      kp, meth_list, chrom,
      meth_r0, meth_r1, "mCpG",
      "#b88690"
    )

    kp <- plot_bigwig(
      kp, h3k4me3_list, chrom,
      h3k4me3_r0, h3k4me3_r1, "H3K4me3",
      "#827A60"
    )

    kp <- plot_bigwig(
      kp, h3k27ac_list, chrom,
      h3k27ac_r0, h3k27ac_r1, "H3K27ac",
      "#848EB5"
    )

    kp <- plot_bigwig(
      kp, h3k27me3_list, chrom,
      h3k27me3_r0, h3k27me3_r1, "H3K27me3",
      "#60474C"
    )
  }
})

dev.off()