# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("karyoploteR", quietly = TRUE, verbose = FALSE)
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
    dplyr::rename(name = motif_id) |>
    dplyr::select(-c(
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

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 chromosome sizes
chm13_chrom_sizes <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0.chrom.sizes"),
  show_col_types = FALSE,
  col_names = c("chr", "size")
) |>
  dplyr::filter(chr != "chrM")

# Add chromosome start (at 0) for each chromosome
chm13_chrom_sizes <- chm13_chrom_sizes |>
  dplyr::mutate(start = 1, .after = chr)

# Load the provided cytobands
chm13_cytobands <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0_cytobands_allchrs.bed"),
  col_names = c("chr", "start", "end", "name", "gieStain"),
  show_col_types = FALSE
)

# Make GRanges objects for both the genome coordinates
# and the cytobands for karyotype plots
chm13_genome_gr <- makeGRangesFromDataFrame(
  chm13_chrom_sizes,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "size",
  keep.extra.columns = TRUE
)

chm13_cytobands_gr <- makeGRangesFromDataFrame(
  chm13_cytobands,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# Load CHM13 segmental duplication annotations
seg_dups <- import(
  fs::path(data, "chm13", "chm13v2.0_SD.bed"),
  format = "BED"
)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

dir_create(c(fs::path(plots, "figure_1")))

pdf(
    fs::path(plots, "figure_1", "karyotype_plot.pdf"),
    width = 8.27,
    height = 11.4
)

pp <- karyoploteR::getDefaultPlotParams(plot.type = 2)
pp$rightmargin <- 0.02
pp$leftmargin <- 0.1
pp$topmargin <- 32
pp$bottommargin <- 32

# Load genome using the T2T-CHM13 genome data
kp <- plotKaryotype(
  genome = chm13_genome_gr,
  cytobands = chm13_cytobands_gr,
  plot.type = 2, # Karyotype with room above for the plot
  cex = 0.9,
  chromosomes = "all",
  plot.params = pp,
  labels.plotter = NULL 
)

kp <- kpAddChromosomeNames(
  kp,
  xoffset = -0.02,
  cex = 0.9
)

# Add base numbers
kp <- kpAddBaseNumbers(
  kp,
  add.units = TRUE,
  cex = 0.5
)

# Dummy density pass to get per-window CTCF counts
kp <- kpPlotDensity(
  kp,
  data = chm13_gr,
  col = NA,
  border = NA,
  window.size = 1e6,
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

# Plot CTCF density per 1Mb
kp <- kpPlotDensity(
  kp,
  chm13_gr,
  window.size = 1e6,
  col = "#c7c8c9",
  data.panel = 1
)

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
    data.panel = 1,
    cex = 0.6,
    tick.pos = chrom_ctcf_max,
    side = 1
  )
}

# Plot segmental duplications
kp <- kpPlotRegions(
  kp,
  GenomicRanges::reduce(seg_dups),
  col = "#D68D00",
  border = NA,
  data.panel = 2,
  r0 = 0.4,
  r1 = 0.7
)

# Plot centromeric and pericentromeric satellites
kp <- kpPlotRegions(
  kp,
  GenomicRanges::reduce(chm13_censat_gr),
  col = "#235EA5",
  border = NA,
  data.panel = 2,
  r0 = 0.7,
  r1 = 1.0
)

dev.off()
