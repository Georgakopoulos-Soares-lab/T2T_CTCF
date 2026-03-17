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
    library("paletteer", quietly = TRUE, verbose = FALSE)
    library("svglite", quietly = TRUE, verbose = FALSE)
    library("zeallot", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores the
# paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

# Path containing FIMO results for CHM13, HG002,
# T2T primates, and human pangenome assemblies
fimo_res <- fs::path(config$fimo_res)

# Path containing metadata for CHM13, HG002,
# T2T primates, and human pangenome assemblies
data <- fs::path(config$data)

plots <- fs::path(config$plots) # Path to store any generated plots

load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
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

chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 chromosome sizes
chm13_chrom_sizes <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0.fa.chrom.sizes"),
  show_col_types = FALSE,
  col_names = c("chr", "size")
)

# Add chromosome start (at 0) for each chromosome
chm13_chrom_sizes <- chm13_chrom_sizes |>
  mutate(start = 0, .after = chr)

# Load the provided cytobands
chm13_cytobands <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0_cytobands_allchrs.bed"),
  col_names = c("chr", "start", "end", "name", "gieStain"),
  show_col_types = FALSE
)

# Make GRanges objects for both the genome
# coordinates and the cytobands for karyotype plots
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

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# A mapping that goes from abbreviations to full names (used for plotting)
censat_mapping <- c(
  "ct" = "Centromeric transitions",
  "censat" = "Other satellites",
  "gsat" = "γSat",
  "bsat" = "βSat",
  "hsat3" = "HSat3",
  "hsat2" = "HSat2",
  "hsat1B" = "HSat1B",
  "hsat1A" = "HSat1A",
  "mon" = "αSat monomeric",
  "dhor" = "αSat divergent",
  "hor" = "αSat inactive HOR",
  "hor_l" = "αSat active HOR",
  "rDNA" = "rDNA array"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- sub("_.*", "", chm13_censat_gr$name)

chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & grepl("L\\)$", chm13_censat_gr$name)
] <- "hor_l"

# Drop unnecessary columns
mcols(chm13_censat_gr)[, c("score", "itemRgb", "thick")] <- NULL

# Karyoplot helper functions
init_karyoplot <- function(
  genome, cytobands, chromosome,
  zoom, plot.type, leftmargin,
  topmargin, bottommargin,
  data1inmargin, ideogramheight,
  cex, tick.dist, minor.tick.dist, digits
) {

  pp <- karyoploteR::getDefaultPlotParams(plot.type = plot.type)
  pp$ideogramheight <- ideogramheight
  pp$leftmargin <- leftmargin
  pp$topmargin <- topmargin
  pp$bottommargin <- bottommargin
  pp$data1inmargin <- data1inmargin

  kp <- karyoploteR::plotKaryotype(
    genome = genome,
    cytobands = cytobands,
    plot.type = plot.type,
    plot.params = pp,
    cex = cex,
    chromosomes = chromosome,
    zoom = zoom
  )

  kp <- karyoploteR::kpAddBaseNumbers(
    kp,
    tick.dist = tick.dist,
    minor.tick.dist = minor.tick.dist,
    add.units = TRUE,
    cex = cex - 0.1,
    digits = digits
  )

  kp
}

plot_density_track <- function(
  kp, data, zoom, r0, r1,
  window.size, data.panel
) {

  kp <- karyoploteR::kpPlotDensity(
    kp, data = data, window.size = window.size,
    r0 = r0, r1 = r0 + 0.01
  )

  computed <- kp$latest.plot$computed.values
  ymax <- max(computed$density[computed$windows %over% zoom])

  kp <- karyoploteR::kpPlotDensity(
    kp, data = data, window.size = window.size,
    ymax = ymax, col = "#c7c8c9", lwd = 1,
    data.panel = data.panel,
    r0 = r0, r1 = r1, clipping = TRUE
  )

  list(kp = kp, ymax = ymax)
}

plot_methylation <- function(
  kp, path, region,
  window.size, r0, r1
) {

  bw <- rtracklayer::BigWigFile(path)

  bins <- rtracklayer::summary(
    bw,
    which = region,
    size = floor(IRanges::width(region) / window.size),
    type = "mean"
  )[[1]]

  kp <- karyoploteR::kpArea(
    kp,
    data = bins,
    y = bins$score,
    ymax = 1,
    col = "#f9d3d3",
    data.panel = 1,
    r0 = r0,
    r1 = r1
  )

  kp
}

add_track_axis <- function(
  kp, ymin, ymax,
  r0, r1, label,
  label_margin, label_pos,
  cex_axis, cex_label
) {

  label_mid <- label_pos %||% (r0 + r1) / 2

  kp <- karyoploteR::kpAxis(
    kp, ymin = ymin, ymax = ymax,
    data.panel = 1, r0 = r0, r1 = r1,
    cex = cex_axis
  )

  kp <- karyoploteR::kpAddLabels(
    kp, labels = label, r0 = label_mid,
    r1 = label_mid, cex = cex_label,
    label.margin = label_margin
  )

  kp
}

add_threshold_lines <- function(
  kp, thresholds, r0, r1,
  data.panel
) {

  for (thr in thresholds) {
    kp <- karyoploteR::kpAbline(
      kp, h = thr$h, ymax = 1,
      r0 = r0, r1 = r1, data.panel = data.panel,
      col = thr$col, lty = 2, lwd = 1
    )
  }

  kp
}

svglite::svglite(
  fs::path(plots, "chr14_acro_zoom1.svg"),
  width = 4.3,
  height = 2.8
)

# Define the region for the chr14 acrocentric short arm
chr14_acro <- GRanges(
  "chr14",
  IRanges(
    1,
    # up until the start of the active alpha satellite
    start(chm13_censat_gr[
      seqnames(chm13_censat_gr) == "chr14" &
        chm13_censat_gr$type == "hor_l"
    ])
  )
)

kp <- init_karyoplot(
  genome = chm13_genome_gr, cytobands = chm13_cytobands_gr,
  chromosome = "chr14", zoom = chr14_acro, plot.type = 2,
  leftmargin = 0.23, topmargin = 0, bottommargin = 0.5,
  data1inmargin = 5, ideogramheight = 15, cex = 0.5,
  tick.dist = 1e6, minor.tick.dist = 5e5, digits = 2
)

window <- 1e5

ctcf_r0 <- 0.00
ctcf_r1 <- 0.45

meth_r0 <- 0.525
meth_r1 <- 1.00

c(kp, ctcf_max) %<-% plot_density_track(
  kp = kp, data = chm13_gr, zoom = chr14_acro,
  r0 = ctcf_r0, r1 = ctcf_r1, window.size = window,
  data.panel = 1
)

kp <- add_track_axis(
  kp = kp, ymin = 0, ymax = ctcf_max,
  r0 = ctcf_r0, r1 = ctcf_r1,
  label = "CTCF density\n(sites/100Kb)",
  label_margin = 0.072, label_pos = NULL,
  cex_axis = 0.45, cex_label = 0.5
)

kp <- plot_methylation(
  kp = kp,
  path = fs::path(
    data, "chm13",
    "chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw"
  ),
  region = chr14_acro, window.size = window,
  r0 = meth_r0, r1 = meth_r1
)

kp <- add_track_axis(
  kp = kp, ymin = 0, ymax = 1,
  r0 = meth_r0, r1 = meth_r1,
  label = "Mean methylation\n(score/100Kb)",
  label_margin = 0.072, label_pos = NULL,
  cex_axis = 0.45, cex_label = 0.5
)

kp <- add_threshold_lines(
  kp = kp,
  thresholds = list(
    list(h = 0.3, col = "#247cd3"),
    list(h = 0.7, col = "#ca2e40")
  ),
  r0 = meth_r0, r1 = meth_r1, data.panel = 1
)

# Define the order of the CenSat types
types <- c(
  "hor_l", "hor", "mon",
  "dhor", "bsat", "gsat",
  "hsat1A", "hsat1B", "hsat2",
  "hsat3", "ct", "censat",
  "rDNA"
)

# Get a color for each type from the paletteer library
colors <- paletteer_d("ggthemes::stata_s2color", n = length(types))
names(colors) <- types

# Go over each CenSat type and plot the regions
for (i in seq_along(types)) {
  region_subset <- chm13_censat_gr[chm13_censat_gr$type == types[i]]
  at <- autotrack(i, length(types), margin = 0.1, r0 = 0.0, r1 = 1.0)

  kp <- kpPlotRegions(
    kp,
    data = region_subset,
    col  = colors[i],
    r0 = at$r0,
    r1 = at$r1,
    border = NA,
    data.panel = 2
  )

  kp <- kpAddLabels(
    kp,
    labels = censat_mapping[types[i]],
    r0 = at$r0,
    r1 = at$r1,
    cex = 0.5,
    data.panel = 2
  )
}
dev.off()

suppressMessages({
  suppressWarnings({
    txdb <- txdbmaker::makeTxDbFromGFF(
      fs::path(data, "chm13", "chm13v2.0_RefSeq_Liftoff_v5.2.gff3"),
      format = "gff3",
      organism = "Homo sapiens",
      taxonomyId = 9606
    )
  })
})

gff3 <- import(fs::path(data, "chm13", "chm13v2.0_RefSeq_Liftoff_v5.2.gff3"))

chr14_rdna_start <- start(
  chm13_censat_gr[
    seqnames(chm13_censat_gr) == "chr14" & chm13_censat_gr$type == "rDNA"
  ]
)

chr14_rdna_end <- end(
  chm13_censat_gr[
    seqnames(chm13_censat_gr) == "chr14" & chm13_censat_gr$type == "rDNA"
  ]
)

chr14_rdna_region <- GRanges("chr14", IRanges(chr14_rdna_start, chr14_rdna_end))

svglite::svglite(
  fs::path(plots, "chr14_acro_zoom2.svg"),
  width = 4.3,
  height = 2.8
)

kp <- init_karyoplot(
  genome = chm13_genome_gr, cytobands = chm13_cytobands_gr,
  chromosome = "chr14", zoom = chr14_rdna_region, plot.type = 1,
  leftmargin = 0.23, topmargin = 0, bottommargin = 10,
  data1inmargin = 2, ideogramheight = 5, cex = 0.5,
  tick.dist = 1e5, minor.tick.dist = 5e4, digits = 2
)

window <- 1e4

rdna_r0 <- 0.0
rdna_r1 <- 0.2

suppressWarnings({
  suppressMessages({
    rdna_genes <- karyoploteR::makeGenesDataFromTxDb(
      txdb = txdb,
      karyoplot = kp,
      plot.transcripts = TRUE,
      plot.transcripts.structure = TRUE
    )

    rna45s2_ids <- as.data.frame(mcols(gff3)) |>
      dplyr::select(c(ID, gene_name)) |>
      filter(gene_name == "RNA45S2") |>
      pull(ID)

    rna45s2_genes <- mcols(rdna_genes$genes)$gene_id %in% rna45s2_ids

    rdna_genes$genes <- rdna_genes$genes[rna45s2_genes]
    rdna_genes$transcripts <- rdna_genes$transcripts[rna45s2_genes]

    kp <- kpPlotGenes(
      kp,
      data = rdna_genes,
      col = "#AD8F8A",
      gene.names = NA,
      r0 = rdna_r0,
      r1 = rdna_r1
    )

    kp <- kpAddLabels(
      kp,
      labels = "rDNA gene copies",
      r0 = 0.07,
      r1 = 0.07,
      cex = 0.5,
      label.margin = 0.05
    )
  })
})

ctcf_r0 <- 0.15
ctcf_r1 <- 0.35

c(kp, ctcf_max) %<-% plot_density_track(
  kp = kp, data = chm13_gr, zoom = chr14_rdna_region,
  r0 = ctcf_r0, r1 = ctcf_r1, window.size = window,
  data.panel = 1
)

kp <- add_track_axis(
  kp = kp, ymin = 0, ymax = ctcf_max,
  r0 = ctcf_r0, r1 = ctcf_r1,
  label = "CTCF density\n(sites/10Kb)",
  label_margin = 0.072, label_pos = NULL,
  cex_axis = 0.45, cex_label = 0.5
)

meth_r0 <- 0.4
meth_r1 <- 0.6

kp <- plot_methylation(
  kp = kp,
  path = fs::path(
    data, "chm13",
    "chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw"
  ),
  region = chr14_rdna_region, window.size = window,
  r0 = meth_r0, r1 = meth_r1
)

kp <- add_track_axis(
  kp = kp, ymin = 0, ymax = 1,
  r0 = meth_r0, r1 = meth_r1,
  label = "Mean methylation\n(score/10Kb)",
  label_margin = 0.072, label_pos = NULL,
  cex_axis = 0.45, cex_label = 0.5
)

kp <- add_threshold_lines(
  kp = kp,
  thresholds = list(
    list(h = 0.3, col = "#247cd3"),
    list(h = 0.7, col = "#ca2e40")
  ),
  r0 = meth_r0, r1 = meth_r1, data.panel = 1
)

dev.off()

pdf(
  fs::path(plots, "chr14_acro_zoom3.pdf"),
  width = 4.3,
  height = 2.8
)

rdna_gene <- GRanges("chr14:2189479-2205433")

kp <- init_karyoplot(
  genome = chm13_genome_gr, cytobands = chm13_cytobands_gr,
  chromosome = "chr14", zoom = rdna_gene, plot.type = 1,
  leftmargin = 0.23, topmargin = 0, bottommargin = 11,
  data1inmargin = 2, ideogramheight = 5, cex = 0.5,
  tick.dist = 1e4, minor.tick.dist = 5e3, digits = 2
)

suppressWarnings({
  rdna_genes <- makeGenesDataFromTxDb(
    txdb,
    karyoplot = kp,
    plot.transcripts = TRUE,
    plot.transcripts.structure = TRUE
  )
})

window <- 1e3

rdna_r0 <- 0.05
rdna_r1 <- 0.3

kp <- kpPlotGenes(
  kp,
  data = rdna_genes,
  col = "#AD8F8A",
  gene.names = NA,
  r0 = rdna_r0,
  r1 = rdna_r1
)

kp <- kpAddLabels(
  kp,
  labels = "rDNA gene copy",
  r0 = rdna_r0,
  r1 = rdna_r1,
  cex = 0.5,
  label.margin = 0.065
)

ctcf_data <- reduce(chm13_gr, ignore.strand = FALSE)

ctcf_cols <- ifelse(
  strand(ctcf_data) == "-", "#c79a46",
  ifelse(strand(ctcf_data) == "+", "#1b9688", "#505152")
)

ctcf_r0 <- 0.3
ctcf_r1 <- 0.35

kp <- kpPlotRegions(
  kp,
  data = ctcf_data,
  col  = ctcf_cols,
  r0 = ctcf_r0,
  r1 = ctcf_r1,
  border = ctcf_cols
)

ctcf_mid <- NULL %||% (ctcf_r0 + ctcf_r1) / 2

kp <- kpAddLabels(
  kp,
  labels = "CTCF-binding\nsites       ",
  r0 = ctcf_mid,
  r1 = ctcf_mid,
  cex = 0.5,
  label.margin = 0.065
)

meth_r0 <- 0.4
meth_r1 <- 0.6

kp <- plot_methylation(
  kp = kp,
  path = fs::path(
    data, "chm13",
    "chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw"
  ),
  region = rdna_gene, window.size = window,
  r0 = meth_r0, r1 = meth_r1
)

kp <- add_track_axis(
  kp = kp, ymin = 0, ymax = 1,
  r0 = meth_r0, r1 = meth_r1,
  label = "Mean methylation\n(score/1Kb)",
  label_margin = 0.072, label_pos = NULL,
  cex_axis = 0.45, cex_label = 0.5
)

kp <- add_threshold_lines(
  kp = kp,
  thresholds = list(
    list(h = 0.3, col = "#247cd3"),
    list(h = 0.7, col = "#ca2e40")
  ),
  r0 = meth_r0, r1 = meth_r1, data.panel = 1
)

dev.off()