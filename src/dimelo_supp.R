# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
    library("EnrichedHeatmap", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
plots_path <- fs::path(res_path, "plots") # Path to store any generated plots

load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "# ", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE,
  ) |>
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

  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data_path, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Load the CTCF DiMeLo-seq aggregated signal
bw_file <- fs::path(res_path, "misc", "CTCF_dimelo.bw")

# Helper function to create the DiMeLo-seq signal plot
plot_regions <- function(
  file,
  regions,
  width,
  extension,
  title_text
) {

  # Convert each CTCF range to its center
  centers <- resize(
    regions,
    width = 1,
    fix = "center"
  )
  # Define a query region (exactly what will be plotted)
  # to load from the DiMeLo-seq BigWig file
  query <- centers + extension

  # Load DiMeLo-seq scores
  bw <- rtracklayer::import(file, which = query)

  # Create EnrichedHeatmap matrix
  mat <- normalizeToMatrix(
    signal = bw,
    target = centers,
    extend = extension,
    w = width,
    value_column = "score",
    background = 0,
    mean_mode = "weighted",
    smooth = FALSE,
    verbose = FALSE
  )

# Create the dataframe we will use for plotting
plot_df <- tibble::tibble(
    bin = seq_len(ncol(mat)),
    mean = colMeans(mat, na.rm = TRUE)
  )

# Extract the indices for where the normalized region starts and ends
# and for where the flanking regions end
n_bins <- ncol(mat)
center_bin <- length(attr(mat, "upstream_index")) + 0.5
breaks = c(1, center_bin, n_bins)

# Create the metaprofile plot and save to a PDF
ggplot(
  plot_df, aes(x = bin)
) +
  geom_line(
    data = plot_df,
    aes(y = mean),
    color = "#3C3C3C", linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = breaks,
    labels = c("-1 Kb", "CTCF site center", "+1 Kb"),
    expand = c(0.05, 0.05)
  ) +
  labs(y = "Mean 6mA (%)", x = NULL, title = title_text) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.background = NULL,
    legend.position = "none"
  )
}

# Load the teucer array annotations
teucer_df <- readr::read_tsv(
  fs::path(data_path, "chm13", "chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"),
  show_col_types = FALSE,
  col_names = c(
   "seqnames", "start", "end",
   "name", "score", "strand",
   "class", "subclass",
   "itemRgb", "blockCount",
   "blockSizes", "blockStarts",
   "RM index", "orig_out"
  )
) |>
  dplyr::select(
    seqnames, start, end, name, strand, class, subclass
  ) |>
  dplyr::filter(
    grepl("teucer", name)
  )

# Convert to GRanges object
teucer_gr <- makeGRangesFromDataFrame(
    teucer_df,
    keep.extra.columns = TRUE
)

# Detect which CTCF sites overlap teucer arrays
teucer_hits <- GenomicRanges::findOverlaps(
  chm13_gr,
  teucer_gr,
  ignore.strand = TRUE
)

# Extract those sites and create the plot
teucer_ctcf_sites <- chm13_gr[queryHits(teucer_hits)]

teucer_plot <- plot_regions(
  file = bw_file,
  regions = teucer_ctcf_sites,
  width = 25,
  extension = 1000,
  title_text = "Genome-wide teucer arrays"
)

ggsave(
    fs::path(plots_path, "dimelo_teucer.pdf"),
    teucer_plot,
    width = 3.6,
    height = 3.6,
    bg = "transparent"
)

# Extract the centromeric transition ranges and
# find which CTCF sites overlap them
ct_annotations <- chm13_censat_gr[chm13_censat_gr$type == "ct"]
ct_hits <- GenomicRanges::findOverlaps(
  chm13_gr,
  ct_annotations,
  ignore.strand = TRUE
)

# Extract centromeric transition overlapping sites and plot
ct_ctcf_sites <- chm13_gr[queryHits(ct_hits)]
ct_plot <- plot_regions(
  file = bw_file,
  regions = ct_ctcf_sites,
  width = 20,
  extension = 1000,
  title_text = "Genome-wide\ncentromeric transition regions"
)

ggsave(
    fs::path(plots_path, "dimelo_ct.pdf"),
    ct_plot,
    width = 3.6,
    height = 3.6,
    bg = "transparent"
)

# Load CHM13 GFF3 file
gff3 <- rtracklayer::import(
  fs::path(data_path, "chm13", "chm13.draft_v2.0.gene_annotation.gff3")
)

# Extract rDNA regions and subset the GFF3 file to contain only
# rRNA genes
rdna_arrays <- chm13_censat_gr[chm13_censat_gr$type == "rdna"]
rdna_gff3 <- subsetByOverlaps(gff3, rdna_arrays, ignore.strand = TRUE)

# Extract RNA45S genes across all chromosomes
rdna_genes <- sort(rdna_gff3[
  rdna_gff3$type == "gene" &
    rdna_gff3$source_gene_common_name %in% c(outer("RNA45S", 1:5, paste0))
])

# Extract rDNA promoter regions and select only chr15 (case study in manuscript)
rdna_promoters <- promoters(rdna_genes)
rdna_promoters <- rdna_promoters[seqnames(rdna_promoters) == "chr15"]

# Find which CTCF sites overlap the chromosome 15 rDNA promoters
rdna_promoter_hits <- GenomicRanges::findOverlaps(
  chm13_gr,
  rdna_promoters,
  ignore.strand = TRUE
)

# Extract those sites and plot
rdna_promoter_ctcf_sites <- chm13_gr[queryHits(rdna_promoter_hits)]
rdna_promoter_plot <- plot_regions(
  file = bw_file,
  regions = rdna_promoter_ctcf_sites,
  width = 20,
  extension = 1000,
  title_text = "Chromosome 15 rDNA promoters"
)

ggsave(
    fs::path(plots_path, "dimelo_rdna_promoters.pdf"),
    rdna_promoter_plot,
    width = 3.6,
    height = 3.6,
    bg = "transparent"
)

# Find which CTCF sites overlap the rDNA arrays in general
rdna_array_hits <- GenomicRanges::findOverlaps(
  chm13_gr,
  rdna_arrays,
  ignore.strand = TRUE
)

# Extract those sites and plot
rdna_array_ctcf_sites <- chm13_gr[
  queryHits(rdna_array_hits)
]

rdna_array_plot <- plot_regions(
  file = bw_file,
  regions = rdna_array_ctcf_sites,
  width = 20,
  extension = 1000,
  title_text = "Genome-wide rDNA arrays"
)

ggsave(
  fs::path(plots_path, "dimelo_rdna.pdf"),
  rdna_array_plot,
  width = 3.6,
  height = 3.6,
  bg = "transparent"
)