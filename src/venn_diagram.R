# Load necessary libraries
suppressPackageStartupMessages({
    suppressWarnings({
        library("rtracklayer", quietly = TRUE, verbose = FALSE)
        library("purrr", quietly = TRUE, verbose = FALSE)
        library("dplyr", quietly = TRUE, verbose = FALSE)
        library("GenomicRanges", quietly = TRUE, verbose = FALSE)
        library("readr", quietly = TRUE, verbose = FALSE)
        library("fs", quietly = TRUE, verbose = FALSE)
        library("yaml", quietly = TRUE, verbose = FALSE)
        library("karyoploteR", quietly = TRUE, verbose = FALSE)
        library("ggvenn", quietly = TRUE, verbose = FALSE)
        library("ggplot2", quietly = TRUE, verbose = FALSE)
        library("patchwork", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res = fs::path(config$res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
plots = fs::path(config$plots) # Path to store any generated plots

motifs <- c("MA0139.2", "MA1929.2", "MA1930.2")

read_beds <- function(prefix, suffix = ".bed") {
  files <- fs::path(res, "misc", paste0(prefix, motifs, suffix))
  grl <- GRangesList(lapply(files, import))
  names(grl) <- motifs
  grl
}

chm13_grl <- read_beds("chm13_")
chm13_liftover_grl <- read_beds("chm13_", ".liftover.hg38.bed")
hg38_grl <- read_beds("hg38_")

# Remove any liftover mappings that are less than the smallest CTCT motif
chm13_liftover_grl <- endoapply(chm13_liftover_grl, function(gr) {
  gr$length <- width(gr)
  gr[width(gr) >= 15]
})

make_venn <- function(motif) {
  # Extract the GRanges object for each motif across hg38, T2T-CHM13, and the liftover annotations
  hg38_motif_gr <- hg38_grl[[motif]]
  chm13_liftover_motif_gr <- chm13_liftover_grl[[motif]]
  chm13_orig_motif_gr <- chm13_grl[[motif]]
  
  # Extract how many regions failed to liftover
  failed_liftover_count <- length(setdiff(unique(chm13_orig_motif_gr$name), unique(chm13_liftover_motif_gr$name)))
  failed_liftover_ids <- paste0("chm13_only_", seq_len(failed_liftover_count))
  
  # Find which hg38 CTCF sites overlap completely (type = within) with any T2T-CHM13 liftover region.
  hits <- findOverlaps(hg38_motif_gr, chm13_liftover_motif_gr, type = "within")
  
  # Assign comparable names between T2T-CHM13 and hg38 for the set operations need in the Venn diagramm
  hg38_ids <- paste0("hg38_", seq_along(hg38_motif_gr))
  chm13_liftover_ids <- paste0("chm13_liftover_", seq_along(chm13_liftover_motif_gr))
  
  # Find which ids are shared between the two (these are the T2T-CHM13 CTCF sites that have a corresponding site in the hg38 genome)
  overlap_hg38_idx <- unique(queryHits(hits))
  overlap_chm13_liftover_idx <- unique(subjectHits(hits))
  shared_ids <- paste0("shared_", seq_along(overlap_hg38_idx))
  
  # Create the sets and assign the unique and shared ids
  hg38_set <- hg38_ids
  hg38_set[overlap_hg38_idx] <- shared_ids
  chm13_liftover_set <- chm13_liftover_ids
  chm13_liftover_set[overlap_chm13_liftover_idx] <- shared_ids

  # In T2T-CHM13 add also the regions that failed to liftover since these are unique to T2T-CHM13
  chm13_liftover_set <- c(chm13_liftover_set, failed_liftover_ids)
  
  # Create the Venn diagramm
  ggvenn(
    list(hg38 = hg38_set, `T2T-CHM13` = chm13_liftover_set),
    fill_color = c("#b64547","#566B94"),
    stroke_size = 0.2,
    set_name_size = 2,
    text_size = 2,
    stroke_alpha = 1,
    fill_alpha = 0.5,
    show_percentage = FALSE
  ) +
    ggtitle(motif) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 7)
    )
}

# Hold all plots in a list
plots_list <- lapply(motifs, function(motif) {
  suppressWarnings(make_venn(motif))
})

# Combine them into a patchwork plot
combined <- wrap_plots(plots_list, ncol = 3)
ggsave(
  fs::path(plots, "all_liftover_venns.pdf"),
  combined, 
  width = 4.5,
  height = 1.5,
  bg = "transparent"
)