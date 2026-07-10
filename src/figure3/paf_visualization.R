# Load required libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("SVbyEye", quietly = TRUE, verbose = FALSE)
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("patchwork", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("circlize", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores the
# paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
plots_path <- fs::path(res_path, "plots") # Path to store any generated plots

# Read the PAF table from, where we have alligned
# the target SD region to the other primate genomes
paf_file <- fs::path(res_path, "misc", "regions_vs_primates.paf")

paf_table <- readPaf(
  paf.file = paf_file,
  include.paf.tags = TRUE
)

# The offset represents the region:chr1_hap1_hsa1:209537154-209792028
# Since we the alignment only used the only the specific segmental duplication
# region, it starts at 0. So we add 209537154 to both start and end to reach the
# genomic position of the target region in the genome
offset <- 209537154

paf_table <- paf_table |>
  dplyr::mutate(
    q.start = q.start + offset,
    q.end = q.end + offset
  )

# Create a name map to convert the PAF species names to
# their scientific notation for plotting
primate_name_map <- c(
  "human" = "Homo sapiens",
  "pon_pyg" = "Pongo pygmaeus",
  "panTro" = "Pan troglodytes"
)

# Create a vector of the species we want to process and show in the plot
species.vec <- c(
  "human",
  "pon_pyg"
)

# This function reads all three CTCF motif FIMO results
# and loads them into one GRanges object
load_fimo_results <- function(path) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(path, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
  ) |>
    select(-c(
      motif_alt_id, score, `p-value`, `q-value`, matched_sequence
    )) # Remove unnecessary columns
  )

  # Create one complete dataframe with the results
  res <- dplyr::bind_rows(ctcf_dfs)

  # Convert the dataframe into a GRanges
  # object for downstream interval manipulation
  ctcf_gr <- GenomicRanges::makeGRangesFromDataFrame(
    res,
    seqnames.field = "sequence_name",
    start.field = "start",
    end.field = "stop",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

}

# Load CTCF FIMO sites for CHM13, Pongo pygmaeus, and Pan troglodytes
chm13_gr <- load_fimo_results(
  fs::path(fimo_res, "chm13")
)

ponpyg_gr <- load_fimo_results(
  fs::path(fimo_res, "primates", "PonPyg")
)

pantro_gr <- load_fimo_results(
  fs::path(fimo_res, "primates", "PanTro")
)

# Define a list with column names used in reading the segmental duplications BED
SD_COLS <- c(
  "chr1", "start1", "end1", "name", "fakeScore", "strand1",
  "thickStart", "thickEnd", "color", "chr2", "start2", "end2",
  "score", "strand2", "max_len", "aln_len", "indel_a", "indel_b",
  "alnB", "matchB", "mismatchB", "transitionsB", "transversions", "fracMatch",
  "fracMatchIndel", "jck", "k2K", "aln_gaps", "uppercaseA", "uppercaseB",
  "uppercaseMatches", "aln_matches", "aln_mismatches", "aln_gap_bases",
  "filter_score", "sat_bases", "unique_id", "original",
  "telo", "peri", "acro", "telo2", "peri2", "acro2"
)

# Define also duplication similarity levels
SIM_LEVELS <- c("<95%", "95-98%", "98-99%", ">=99%")

# Τhis is a helper function that reads the PAF table, filters to
# a specific species, and converts to a GRanges object
convert_paf_table_to_gr <- function(
  paf_table,
  paf_species_name,
  prefix
) {
  name_col <- paste0(prefix, "name")
  chr_col <- paste0(prefix, "chr")
  start_col <- paste0(prefix, "start")
  end_col <- paste0(prefix, "end")

  out <- paf_table |>
    dplyr::filter(
      grepl(paf_species_name, .data[[name_col]], ignore.case = TRUE)
    ) |>
    dplyr::mutate(
      !!chr_col := sub(":[0-9]+-[0-9]+$", "", sub("^[^#]*#", "", .data[[name_col]]))
    )

  GenomicRanges::makeGRangesFromDataFrame(
    out,
    seqnames.field = chr_col,
    start.field = start_col,
    end.field = end_col,
    keep.extra.columns = TRUE
  )
}

convert_sd_to_gr <- function(
  sd_full_df,
  paf_table,
  paf_species_name,
  prefix
) {
  # Convert SD regions to a GRanges object
  sd_gr <- GRanges(
    seqnames = sd_full_df$seqnames,
    ranges   = IRanges(start = sd_full_df$start, end = sd_full_df$end),
    strand   = sd_full_df$strand,
    similarity = sd_full_df$similarity,
    sim_class = sd_full_df$sim_class
  )

  paf_gr <- convert_paf_table_to_gr(paf_table, paf_species_name, prefix)

  # Keep only SD coordinates that overlap with the PAF regions of interest
  filtered_sd_gr <- IRanges::subsetByOverlaps(
    sd_gr,
    paf_gr,
    ignore.strand = TRUE
  )

  # Create a final GRanges object that uses as seqnames the species name
  # so that it aligns with the PAF table used in SVbyEye
  final_sd_gr <- GRanges(
    seqnames = primate_name_map[[paf_species_name]],
    ranges = ranges(filtered_sd_gr),
    strand = strand(filtered_sd_gr)
  )

  # Update metadata columns and return the final object 
  mcols(final_sd_gr) <- mcols(filtered_sd_gr)

  final_sd_gr
}

# Define a function that reads the segmental duplication BED
# filters to keep only the columns needed, add similarity categories
# converts to GRanges object for ploting 
read_sd <- function(
  filepath,
  paf_table,
  paf_species_name,
  prefix
) {
  sd_df <- readr::read_tsv(
    filepath,
    col_names = SD_COLS,
    col_select = c(chr1, start1, end1, strand1, fracMatch),
    show_col_types = FALSE
  ) |>
  dplyr::rename(
    seqnames = chr1,
    start = start1,
    end = end1,
    strand = strand1,
    similarity = fracMatch
  ) |>
  dplyr::mutate(
    sim_class = factor(
      dplyr::case_when(
        similarity >= 0.99 ~ ">=99%",
        similarity >= 0.98 ~ "98-99%",
        similarity >= 0.95 ~ "95-98%",
        .default = "<95%"
      ),
      levels = SIM_LEVELS,
      ordered = TRUE
    )
  )

  convert_sd_to_gr(
    sd_df,
    paf_table,
    paf_species_name,
    prefix
  )
}

# Load CHM13 SD regions
chm13_sd_gr <- read_sd(
  fs::path(data_path, "chm13", "chm13v2.0_SD.full.bed"),
  paf_table,
  paf_species_name = "human",
  prefix = "t."
)

ponpyg_sd_gr <- read_sd(
  fs::path(data_path, "primates", "PonPyg", "mPonPyg2_v2.0.SD_v1.0.bed"),
  paf_table,
  paf_species_name = "pon_pyg",
  prefix = "t."
)

pantro_sd_gr <- read_sd(
  fs::path(data_path, "primates", "PanTro", "mPanTro3_v2.0.SD_v1.0.bed"),
  paf_table,
  paf_species_name = "panTro",
  prefix = "q."
)

# This is a helper function that will filter CTCF sites to those only
# within the alignment region
filter_ctcf <- function(
  ctcf_sites,
  paf_table,
  paf_species_name,
  prefix
) {
  
  paf_gr <- convert_paf_table_to_gr(
    paf_table,
    paf_species_name,
    prefix
  )

  filtered_ctcf <- IRanges::subsetByOverlaps(
    ctcf_sites,
    paf_gr,
    ignore.strand = TRUE
  )

  # Create a final GRanges object that uses as seqnames the species name
  # so that it aligns with the PAF table used in SVbyEye
  final_ctcf_gr <- GRanges(
    seqnames = primate_name_map[[paf_species_name]],
    ranges = ranges(filtered_ctcf),
    strand = strand(filtered_ctcf)
  )

  # Update metadata columns and return the final object 
  mcols(final_ctcf_gr) <- mcols(filtered_ctcf)

  final_ctcf_gr

}

# Filter the CTCF sites GRanges object to just CTCF sites within
# the regions present in the PAF file for all three species
chm13_gr <- filter_ctcf(
  chm13_gr,
  paf_table,
  paf_species_name = "human",
  prefix = "t."
)

pantro_gr <- filter_ctcf(
  pantro_gr,
  paf_table,
  paf_species_name = "panTro",
  prefix = "q."
)

ponpyg_gr <- filter_ctcf(
  ponpyg_gr,
  paf_table,
  paf_species_name = "pon_pyg",
  prefix = "t."
)

# Define color palettes for both segmental duplications 
# and for CTCF motifs
sd_palette <- c(
  "<95%"   = "#66006F",
  "95-98%" = "#6D6D6D",
  "98-99%" = "#FCFF00",
  ">=99%"  = "#F29200"
)

motif_palette <- c(
  "MA0139.2" = "#A3ABC0",
  "MA1929.2" = "#424D6A",
  "MA1930.2" = "#0E1933"
)

# Here, we define 2 lists with the annotation objects that will plotted for the query region
# SDs colored by similarity and CTCF sites colored by ID
query_tracks <- list(
  list(gr = pantro_sd_gr, fill.by = "sim_class", palette = sd_palette),
  list(gr = pantro_gr, fill.by = "motif_id", palette = motif_palette)
)

# Here, we define 2 lists with the annotation objects that will plotted for the target regions
# Then we build the 
target_tracks <- list(
  human = list(sd = chm13_sd_gr, ctcf = chm13_gr),
  pon_pyg = list(sd = ponpyg_sd_gr, ctcf = ponpyg_gr)
)

build_specs <- function(species) {
  t <- target_tracks[[species]]
  c(
    list(
      list(gr = t$sd, fill.by = "sim_class", palette = sd_palette, space = "target"),
      list(gr = t$ctcf, fill.by = "motif_id", palette = motif_palette, space = "target")
    ),
    lapply(query_tracks, \(s) c(s, list(space = "query")))
  )
}

# Helper function that modifies the axis so that it is thinner for plotting 
thin_axis <- function(obj, n = 3) {
  lab <- obj$labels
  
  k <- length(lab)
  
  idx  <- unique(round(seq(1, k, length.out = n)))
  orig <- obj$breaks

  obj$breaks <- function(limits, ...) {
    b <- if (is.function(orig)) orig(limits, ...) else orig
    if (length(b) == k) b[idx] else b
  }
  obj$labels <- lab[idx]
  invisible(obj)
}

# Helper function that adds annotation tracks to either target or query species
# This also adjusts the tracks size after plotting since SVbyEye does not have
# size modification flags
add_track <- function(p, spec, blk) {
  n_before <- length(p$layers)

  p <- addAnnotation(
    ggplot.obj = p,
    annot.gr = spec$gr,
    fill.by = spec$fill.by,
    color.palette = spec$palette,
    coordinate.space = spec$space,
    annotation.label = NULL,
    offset.annotation = FALSE
  )

  new_layers <- seq.int(n_before + 1L, length(p$layers))
  for (i in new_layers) {
    p$layers[[i]]$geom_params$arrowhead_height <- blk
    p$layers[[i]]$geom_params$rect_height <- blk
  }

  p
}

# Helper function that sets the same theme across plots
# Theme is standard across plots in the manuscript
panel_theme <- function() {
  list(
    theme_classic(),
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 8, face = "italic", color = "black"),
      axis.text.x = element_text(size = 6),
      axis.title.x = element_text(size = 8),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.5),
      axis.line.x = element_line(linewidth = 0.5),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0.15, unit = "cm")
    )
  )
}

# This is the main function that builds the miro plot for each species
# alongside the added tracks
build_panel <- function(species, paf_table, blk = grid::unit(1, "mm")) {

  # We filter the PAF table to keep only alignments with MAPQ > 10
  # and modify the query and target names to represent species
  d <- paf_table |>
    filter(mapq > 10, grepl(species, t.name)) |>
    separate_wider_delim(t.name, "#", names = c("t.name", "t.chr")) |>
    separate_wider_delim(q.name, "#", names = c("q.name", "q.chr")) |>
    mutate(
      q.name = primate_name_map[q.name],
      t.name = primate_name_map[t.name]
    )

  # We create the miro plot with alignment direction shaded
  p <- plotMiro(
    d,
    color.by = "direction",
    offset.alignments = TRUE,
    genomic.scale = "kbp",
    outline.alignments = TRUE,
    add.alignment.arrows = TRUE,
    color.palette = c("+" = "#1b9688", "-" = "#c79a46")
  ) + panel_theme()

  # Modify the aesthetics of the plot
  p$layers[[2]]$aes_params$colour_ggnewscale_1 <- "black"
  p$layers[[2]]$aes_params$linewidth <- 0.15
  p$layers[[3]]$geom_params$arrow_body_height <- blk
  p$layers[[3]]$geom_params$arrowhead_height  <- blk

  sx <- p$scales$get_scales("x")
  thin_axis(sx, 4)
  if (inherits(sx$secondary.axis, "AxisSecondary")) thin_axis(sx$secondary.axis, 4)

  # We add each target/query annotation
  p_acc <- p

  for (spec in build_specs(species)) {
    p_acc <- add_track(p_acc, spec, blk)
  }

  p_acc
}

# Create the plot for each species
plots <- lapply(species.vec, build_panel, paf_table = paf_table)

# Wrap them into one plot with the shared x-axis title
combined_plot <- wrap_plots(
  plots, 
  ncol = 2,
  axis_titles = "collect"
) & 
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Save plot
ggsave(
  fs::path(plots_path, "primate_paf_alignment.pdf"),
  combined_plot,
  bg = "transparent",
  height = 2.4,
  width = 4.9139
)