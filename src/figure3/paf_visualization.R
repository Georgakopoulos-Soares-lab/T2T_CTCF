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

# Create a name map to convert the PAF species names to
# their scientific notation for plotting
primate_name_map <- c(
  "human" = "Homo sapiens",
  "pan_pan" = "Pan paniscus",
  "pon_abe" = "Pongo abelii",
  "pon_pyg" = "Pongo pygmaeus",
  "panTro" = "Pan troglodytes"
)

# Create a vector of the species we want to process and show in the plot
species.vec <- c(
  "human"
  # "pon_pyg"
)

# Create an empty list that will st
plots <- vector("list", length(species.vec))
names(plots) <- species.vec

read_sd <- function(
  filepath
) {
  suppressMessages({
    sd_full <- read_tsv(
      fs::path(
        filepath
      ),
      show_col_types = FALSE,
      col_names = c(
        "chr1", "start1", "end1", "name", "fakeScore", "strand1",
        "start1", "end1", "color", "chr2", "start2", "end2",
        "score", "strand2", "max_len", "aln_len", "indel_a", "indel_b",
        "alnB", "matchB", "mismatchB", "transitionsB", "transversions", "fracMatch",
        "fracMatchIndel", "jck", "k2K", "aln_gaps", "uppercaseA", "uppercaseB",
        "uppercaseMatches", "aln_matches", "aln_mismatches", "aln_gap_bases", "filter_score",
        "sat_bases", "unique_id", "original", "telo", "peri", "acro", "telo2", "peri2", "acro2"
      )
    ) |>
      dplyr::select(
        `chr1`, `start1...2`, `end1...3`,
        strand1, fracMatch
      ) |>
      dplyr::rename(
        seqnames = `chr1`,
        start = `start1...2`,
        end = `end1...3`,
        strand = strand1,
        similarity = fracMatch
      )|>
      dplyr::mutate(
        sim_class = dplyr::case_when(
          similarity >= 0.99 ~ ">=99%",
          similarity >= 0.98 ~ "98-99%",
          similarity >= 0.95 ~ "95-98%",
          .default = "<95%"
        )
      )
  })

  sd_full
}

# Load CHM13 SD regions
chm13_sd_full <- read_sd(
  fs::path(data_path, "chm13", "chm13v2.0_SD.full.bed")
)

# Define a list with colors for each segmental duplication's
# similarity metric
col_fun <- c(
  "<95%"   = "#66006F",
  "95-98%" = "#6D6D6D",
  "98-99%" = "#FCFF00",
  ">=99%"  = "#F29200"
)

convert_sd_to_gr <- function(
  sd_full_df,
  paf_gr,

) {
  # Convert SD regions to a GRanges object
  sd_gr <- GRanges(
    seqnames = sd_full_df$seqnames,
    ranges   = IRanges(start = sd_full_df$start, end = sd_full_df$end),
    strand   = sd_full_df$strand,
    similarity = sd_full_df$similarity,
    sim_class = sd_full_df$sim_class
  )

  # Keep only SD coordinates that overlap with the PAF regions of interest
  filtered_sd_gr <- GenomicRanges::subsetByOverlaps(
    sd_gr,
    paf_gr,
    ignore.strand = TRUE
  )

  final_sd_gr <- GRanges(
    seqnames = "Homo sapiens",
    ranges = ranges(filtered_sd_gr),
    strand = strand(filtered_sd_gr)
  )

  mcols(final_sd_gr) <- mcols(filtered_sd_gr)

  final_sd_gr
}



paf_gr <- makeGRangesFromDataFrame(
  paf_table |>
    dplyr::filter(grepl("human", t.name)) |>
    separate_wider_delim(
          t.name,
          delim = "#",
          names = c("t.name", "t.chr")
        ),
  seqnames.field = "t.chr",
  start.field = "t.start",
  end.field = "t.end",
  keep.extra.columns = TRUE
)

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

for (species in species.vec) {
    filtered_paf_table <- paf_table |>
      filter(
        mapq >= 30,
        grepl(species, t.name)
      ) 
    
    filtered_paf_table <- filtered_paf_table |>
      separate_wider_delim(
        t.name,
        delim = "#",
        names = c("t.name", "t.chr")
      ) |>
      separate_wider_delim(
        q.name,
        delim = "#",
        names = c("q.name", "q.chr")
      )

    filtered_paf_table$q.name <- primate_name_map[filtered_paf_table$q.name]
    filtered_paf_table$t.name <- primate_name_map[filtered_paf_table$t.name]
    
    if (species == "pon_abe") {
      filtered_paf_table <- filtered_paf_table |>
        dplyr::filter(t.chr == "chr1_hap1_hsa1")
    }

    p <- plotMiro(
      filtered_paf_table,
      color.by = "direction",
      offset.alignments = TRUE,
      genomic.scale = "kbp",
      outline.alignments = TRUE,
      add.alignment.arrows = TRUE,
      color.palette = c(
        "+" = "#1b9688",
        "-" = "#c79a46"
      )
    ) +
      theme_classic() + 
      theme(
        legend.position = "none",
        axis.text.y = element_text(size = 8, face = "italic", color = "black"),
        axis.text.x = element_text(size = 6),
        axis.title.x= element_text(size = 8),
        axis.line.y = element_blank(),
        axis.title.y= element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.line.x  = element_line(linewidth = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0.15, unit = "cm")
      )

    p$layers[[2]]$aes_params$colour_ggnewscale_1 <- "black"
    p$layers[[2]]$aes_params$linewidth <- 0.15
    
    blk <- grid::unit(1.5, "mm")
    
    p$layers[[3]]$geom_params$arrow_body_height <- blk
    p$layers[[3]]$geom_params$arrowhead_height <- blk

    sx <- p$scales$get_scales("x")
    
    thin_axis(sx, 4)
    
    if (inherits(sx$secondary.axis, "AxisSecondary")) {
      thin_axis(sx$secondary.axis, 4)
    }

    p <- addAnnotation(
        ggplot.obj = p,
        annot.gr = chm13_sd_gr,
        fill.by = "sim_class",
        color.palette = col_fun,
        coordinate.space = "target",
        annotation.label = "SDs"
    )

    plots[[species]] <- p
}

combined_plot <- wrap_plots(
  plots, 
  ncol = 2,
  axis_titles = "collect"
) & 
  theme(
    plot.background = element_blank(),
    panel.background = element_blank()
  )

combined_plot

ggsave(
  fs::path(plots_path, "test.pdf"),
  combined_plot,
  bg = "transparent",
  height = 2.4,
  width = 4.9139
)