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
    library("svglite", quietly = TRUE, verbose = FALSE)
    library("zeallot", quietly = TRUE, verbose = FALSE)
  })
})

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

  # Return the FIMO results
  fimo_gr
}

# Helper function which initializes the karyoplot
init_karyoplot <- function(
  genome, chromosome, censat_track,
  zoom, plot.type, leftmargin,
  topmargin, bottommargin,
  data1inmargin, ideogramheight,
  cex, tick_dist, minor_tick_dist, digits
) {

  # Modify the plot parameters for better visulization
  pp <- karyoploteR::getDefaultPlotParams(plot.type = plot.type)
  pp$ideogramheight <- ideogramheight
  pp$leftmargin <- leftmargin
  pp$topmargin <- topmargin
  pp$bottommargin <- bottommargin
  pp$data1inmargin <- data1inmargin

  # Plot karyotype
  kp <- karyoploteR::plotKaryotype(
    genome = genome,
    ideogram.plotter = NULL,
    plot.type = plot.type,
    plot.params = pp,
    cex = cex,
    chromosomes = chromosome,
    zoom = zoom
  )

  # Plot CenSat regions instead of cytobands
  kp <- karyoploteR::kpPlotRegions(
    kp,
    censat_track,
    r0 = 0, r1 = 1,
    col = censat_track$itemRgb,
    border = NA,
    data.panel = "ideogram",
    clipping = TRUE
  )

  # Add base numbers
  kp <- karyoploteR::kpAddBaseNumbers(
    kp,
    tick.dist = tick_dist,
    minor.tick.dist = minor_tick_dist,
    add.units = TRUE,
    cex = cex - 0.1,
    digits = digits
  )

  # Return the plot
  kp
}

# Helper function to plot CTCF density
plot_density_track <- function(
  kp, data, zoom, r0, r1,
  window.size, data.panel
) {

  # Invisible pass to recover per-window densities for ymax
  kp <- karyoploteR::kpPlotDensity(
    kp, data = data, window.size = window.size,
    border = NA, col = NA,
    r0 = r0, r1 = r0 + 0.01
  )

  # Find the density values in the defined window
  computed <- kp$latest.plot$computed.values

  # Extract the maximum density value
  ymax <- max(computed$density[computed$windows %over% zoom])

  # Plot density
  kp <- karyoploteR::kpPlotDensity(
    kp, data = data, window.size = window.size,
    ymax = ymax, border = NA, col = "#3C3C3C",
    lwd = 1, data.panel = data.panel,
    r0 = r0, r1 = r1, clipping = TRUE
  )

  # Return both the plot and the max density for axis placement
  list(kp = kp, ymax = ymax)
}

# Helper function to plot BigWig signals over a specific region
plot_bigwig <- function(kp, path, region, window.size, r0, r1, color) {

  # Create a pointer to the BigWig file
  bw <- rtracklayer::BigWigFile(path)

  # Calculate BigWig score across each bin
  bins <- rtracklayer::summary(
    bw, which = region,
    size = floor(IRanges::width(region) / window.size),
    type = "mean"
  )[[1]]

  # For bins missing signal convert to 0
  bins$score[is.na(bins$score)] <- 0

  # Extract the max score over the window
  ymax <- max(bins$score)

  # Plot
  kp <- karyoploteR::kpArea(
    kp, data = bins, y = bins$score,
    ymax = ymax, col = color, border = color,
    data.panel = 1, r0 = r0, r1 = r1
  )

  # Return both the plot and the max score
  list(kp = kp, ymax = ymax)
}

# Helper function to plot axis
add_track_axis <- function(
  kp, ymin, ymax,
  r0, r1, label,
  label_margin, label_pos,
  cex_axis, cex_label
) {

  # Calculate the axis mid for the label to be plotted
  label_mid <- label_pos %||% (r0 + r1) / 2

  # Plot axis
  kp <- karyoploteR::kpAxis(
    kp, ymin = ymin, ymax = ymax,
    data.panel = 1, r0 = r0, r1 = r1,
    cex = cex_axis, tick.pos = ymax
  )

  # Add labels
  kp <- karyoploteR::kpAddLabels(
    kp, labels = label, r0 = label_mid,
    r1 = label_mid, cex = cex_label,
    label.margin = label_margin
  )

  # Return the updated plot
  kp
}

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 chromosome sizes
chm13_chrom_sizes <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0.chrom.sizes"),
  show_col_types = FALSE,
  col_names = c("chr", "size")
)

# Add chromosome start (at 0) for each chromosome
chm13_chrom_sizes <- chm13_chrom_sizes |>
  mutate(start = 1, .after = chr)

# Make GRanges objects for the genome coordinates
chm13_genome_gr <- makeGRangesFromDataFrame(
  chm13_chrom_sizes,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "size",
  keep.extra.columns = TRUE
)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define alpha satellite HOR as the region that has the associated color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Extract the coordinates for chr15 acrocentric short arm
chr15_acro <- GRanges(
  "chr15",
  IRanges(
    1,
    # including the active alpha satellite
    end(chm13_censat_gr[
      seqnames(chm13_censat_gr) == "chr15" &
        chm13_censat_gr$type == "hor_l"
    ])
  )
)

# CenSat annotations restricted to the acrocentric short arm
chr15_acro_censat <- subsetByOverlaps(
  chm13_censat_gr,
  chr15_acro,
  ignore.strand = TRUE
)

# Define a list with tracks to plot
tracks <- tibble::tribble(
  ~file, ~color, ~label,
  #FiberSeq accessibility
  "all.percent.accessible.bw", "#608266", "% accessible",
  #Nanopore methylation
  "chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw", "#b88690", "mCpG",
  #DiMeLo-seq H3K4me3
  "h3k4me3_dimelo.bw", "#827A60", "H3K4me3",
  #DiMeLo-seq H3K27ac
  "h3k27ac_dimelo.bw", "#848EB5", "H3K27ac",
  #DiMeLo-seq H3K27me3
  "h3k27me3_dimelo.bw", "#60474C", "H3K27me3"
)

# Create plot for zoom 1
# svglite::svglite(
pdf(
  fs::path(plots, "chr15_acro_zoom1.pdf"),
  width = 4.3,
  height = 2
)

# Initialize karyoplot
kp <- init_karyoplot(
  genome = chm13_genome_gr,
  chromosome = "chr15",
  censat_track = chr15_acro_censat,
  zoom = chr15_acro,
  plot.type = 1,
  leftmargin = 0.2,
  topmargin = 0,
  bottommargin = 15,
  data1inmargin = 2,
  ideogramheight = 5,
  cex = 0.5,
  tick_dist = 5e6,
  minor_tick_dist = 1e6,
  digits = 2
)

# In zoom1, we average over 100Kb windows
zoom1_window <- 1e5

# Define gap between tracks and track positioning
# so that they cover the whole plotting region
track_gap <- 0.02
n_tracks <- nrow(tracks) + 1 # +1 for CTCF
track_height <- (1 - (n_tracks - 1) * track_gap) / n_tracks

# Define CTCF track regions
ctcf_r0 <- 0
ctcf_r1 <- track_height

# Define regions to plot the other tracks
tracks$r0 <- ctcf_r1 + track_gap +
  (seq_len(nrow(tracks)) - 1) * (track_height + track_gap)

tracks$r1 <- tracks$r0 + track_height

# Plot CTCF
c(kp, ctcf_max) %<-% plot_density_track(
  kp = kp,
  data = chm13_gr,
  zoom = chr15_acro,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  window.size = zoom1_window,
  data.panel = 1
)

# Add CTCF axis
kp <- add_track_axis(
  kp = kp,
  ymin = 0,
  ymax = ctcf_max,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  label = "CTCF density\n(sites/100Kb)",
  label_margin = 0.072,
  label_pos = NULL,
  cex_axis = 0.45,
  cex_label = 0.5
)

# Iterate over tracks and plot
for (i in seq_len(nrow(tracks))) {
  spec <- tracks[i, ]

  c(kp, track_max) %<-% plot_bigwig(
    kp = kp,
    path = fs::path(data, "chm13", spec$file),
    region = chr15_acro,
    window.size = zoom1_window,
    r0 = spec$r0,
    r1 = spec$r1,
    color = spec$color
  )

  kp <- add_track_axis(
    kp = kp,
    ymin = 0,
    ymax = track_max,
    r0 = spec$r0,
    r1 = spec$r1,
    label = spec$label,
    label_margin = 0.072,
    label_pos = NULL,
    cex_axis = 0.45,
    cex_label = 0.5
  )
}

dev.off()


# Create TxDB object using the Liftoff gene annotations for CHM13
suppressMessages({
  suppressWarnings({
    txdb <- txdbmaker::makeTxDbFromGFF(
      fs::path(data, "chm13", "chm13.draft_v2.0.gene_annotation.gff3"),
      format = "gff3",
      organism = "Homo sapiens",
      taxonomyId = 9606
    )
  })
})

# Load also the associated GFF3 file
gff3 <- import(fs::path(data, "chm13", "chm13.draft_v2.0.gene_annotation.gff3"))

# Extract coordinates for the rDNA region in chr15
chr15_rdna_start <- start(
  chm13_censat_gr[
    seqnames(chm13_censat_gr) == "chr15" & chm13_censat_gr$type == "rdna"
  ]
)

chr15_rdna_end <- end(
  chm13_censat_gr[
    seqnames(chm13_censat_gr) == "chr15" & chm13_censat_gr$type == "rdna"
  ]
)

chr15_rdna_region <- GRanges("chr15", IRanges(chr15_rdna_start, chr15_rdna_end))

# Select only genes and transcripts present within the rDNA region
suppressWarnings({
  suppressMessages({
    rdna_genes <- karyoploteR::makeGenesDataFromTxDb(
      txdb = txdb,
      karyoplot = kp,
      plot.transcripts = TRUE,
      plot.transcripts.structure = TRUE
    )

    rna45s3_ids <- as.data.frame(mcols(gff3)) |>
      dplyr::select(c(ID, gene_name)) |>
      filter(gene_name == "RNA45S3") |>
      pull(ID)

    rna45s3_genes <- mcols(rdna_genes$genes)$gene_id %in% rna45s3_ids

    rdna_genes$genes <- rdna_genes$genes[rna45s3_genes]
    rdna_genes$transcripts <- rdna_genes$transcripts[rna45s3_genes]
  })
})

# Create plot for zoom2
# svglite::svglite(
pdf(
  fs::path(plots, "chr15_acro_zoom2.pdf"),
  width = 4.3,
  height = 2
)

# Initialize plot
kp <- init_karyoplot(
  genome = chm13_genome_gr,
  chromosome = "chr15",
  censat_track = chr15_acro_censat,
  zoom = chr15_rdna_region,
  plot.type = 1,
  leftmargin = 0.23,
  topmargin = 0,
  bottommargin = 15,
  data1inmargin = 2,
  ideogramheight = 5,
  cex = 0.5,
  tick_dist = 1e6,
  minor_tick_dist = 5e5,
  digits = 2
)

# In zoom2 the windows are 10Kb
zoom2_window <- 1e4
rdna_r0 <- 0.0
rdna_r1 <- 0.15

# Plot rDNA genes
suppressWarnings({
  suppressMessages({
    kp <- kpPlotGenes(
      kp,
      data = rdna_genes,
      col = "#AD8F8A",
      gene.names = NA,
      r0 = rdna_r0,
      r1 = rdna_r1
    )
  })
})

# Add rDNA label
kp <- kpAddLabels(
  kp,
  labels = "rDNA gene copies",
  r0 = 0.07,
  r1 = 0.07,
  cex = 0.5,
  label.margin = 0.05
)

# Here, we account for the rDNA track in the available space
available <- 1 - rdna_r1 - track_gap # space above rDNA, minus separator
track_height <- (available - (n_tracks - 1) * track_gap) / n_tracks

# Define regions to plot CTCF and the signal tracks
ctcf_r0 <- 0.09 + track_gap
ctcf_r1 <- ctcf_r0 + track_height

tracks$r0 <- ctcf_r1 + track_gap +
  (seq_len(nrow(tracks)) - 1) * (track_height + track_gap)

tracks$r1 <- tracks$r0 + track_height

# Plot CTCF track and axis
c(kp, ctcf_max) %<-% plot_density_track(
  kp = kp,
  data = chm13_gr,
  zoom = chr15_rdna_region,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  window.size = zoom2_window,
  data.panel = 1
)

kp <- add_track_axis(
  kp = kp,
  ymin = 0,
  ymax = ctcf_max,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  label = "CTCF density\n(sites/10Kb)",
  label_margin = 0.072,
  label_pos = NULL,
  cex_axis = 0.45,
  cex_label = 0.5
)

# Plot BigWig signals and axis
for (i in seq_len(nrow(tracks))) {
  spec <- tracks[i, ]

  c(kp, track_max) %<-% plot_bigwig(
    kp = kp,
    path = fs::path(data, "chm13", spec$file),
    region = chr15_rdna_region,
    window.size = zoom2_window,
    r0 = spec$r0,
    r1 = spec$r1,
    color = spec$color
  )

  kp <- add_track_axis(
    kp = kp,
    ymin = 0,
    ymax = track_max,
    r0 = spec$r0,
    r1 = spec$r1,
    label = spec$label,
    label_margin = 0.072,
    label_pos = NULL,
    cex_axis = 0.45,
    cex_label = 0.5
  )
}

dev.off()

# Reload the rDNA TxDB
suppressWarnings({
  suppressMessages({
    rdna_genes <- karyoploteR::makeGenesDataFromTxDb(
      txdb = txdb,
      karyoplot = kp,
      plot.transcripts = TRUE,
      plot.transcripts.structure = TRUE
    )
  })
})

# Finally, we plot the zoom 3 region which is a single rDNA gene
pdf(
  fs::path(plots, "chr15_acro_zoom3.pdf"),
  width = 4.3,
  height = 2
)

# Coordinates for a single rDNA gene and flanking region
rdna_gene <- GRanges("chr15:2649951-2723289")

# Initialize plot
kp <- init_karyoplot(
  genome = chm13_genome_gr,
  chromosome = "chr15",
  censat_track = chr15_acro_censat,
  zoom = rdna_gene,
  plot.type = 1,
  leftmargin = 0.23,
  topmargin = 0,
  bottommargin = 15,
  data1inmargin = 1,
  ideogramheight = 5,
  cex = 0.5,
  tick_dist = 1e4,
  minor_tick_dist = 5e3,
  digits = 2
)

# In zoom3 the window is 100bp
zoom3_window <- 100
rdna_r0 <- 0.01
rdna_r1 <- 0.2

# Plot rDNA genes
suppressWarnings({
  suppressMessages({
    kp <- kpPlotGenes(
      kp,
      data = rdna_genes,
      col = "#AD8F8A",
      gene.names = NA,
      r0 = rdna_r0,
      r1 = rdna_r1
    )
  })
})

# Add rDNA gene labels
kp <- kpAddLabels(
  kp,
  labels = "rDNA gene copies",
  r0 = 0.07,
  r1 = 0.07,
  cex = 0.5,
  label.margin = 0.05
)

available <- 1 - rdna_r1 - track_gap # space above rDNA, minus separator
track_height <- (available - (n_tracks - 1) * track_gap) / n_tracks

ctcf_r0 <- 0.15 + track_gap
ctcf_r1 <- ctcf_r0 + track_height

tracks$r0 <- ctcf_r1 + track_gap +
  (seq_len(nrow(tracks)) - 1) * (track_height + track_gap)
tracks$r1 <- tracks$r0 + track_height

# Merge overlapping CTCF sites for better plotting
ctcf_data <- reduce(chm13_gr, ignore.strand = FALSE)

# Define a color for + and - strand sites
ctcf_cols <- ifelse(
  strand(ctcf_data) == "-", "#c79a46",
  ifelse(strand(ctcf_data) == "+", "#1b9688", "#505152")
)

# Instead of density plot individual regions
kp <- kpPlotRegions(
  kp,
  data = ctcf_data,
  col  = ctcf_cols,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  border = ctcf_cols
)

# Find the center of the axis for label placement
ctcf_mid <- NULL %||% (ctcf_r0 + ctcf_r1) / 2

kp <- kpAddLabels(
  kp,
  labels = "CTCF-binding\nsites       ",
  r0 = ctcf_mid,
  r1 = ctcf_mid,
  cex = 0.5,
  label.margin = 0.065
)

for (i in seq_len(nrow(tracks))) {
  spec <- tracks[i, ]

  c(kp, track_max) %<-% plot_bigwig(
    kp = kp,
    path = fs::path(data, "chm13", spec$file),
    region = rdna_gene,
    window.size = zoom3_window,
    r0 = spec$r0,
    r1 = spec$r1,
    color = spec$color
  )

  kp <- add_track_axis(
    kp = kp,
    ymin = 0,
    ymax = track_max,
    r0 = spec$r0,
    r1 = spec$r1,
    label = spec$label,
    label_margin = 0.072,
    label_pos = NULL,
    cex_axis = 0.45,
    cex_label = 0.5
  )
}

dev.off()