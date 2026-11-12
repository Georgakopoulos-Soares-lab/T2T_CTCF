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
    library("ggsci", quietly = TRUE, verbose = FALSE)
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

# Define a color palette to be used throughout the plots
motif_palette <- c(
  "MA0139.2" = "#A3ABC0",
  "MA1929.2" = "#424D6A",
  "MA1930.2" = "#0E1933"
)

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

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Add a column specifing the region's length
chm13_censat_gr$length <- width(chm13_censat_gr)

# Loop over the CTCF motifs and find how many overlap
# with the CenSat annotations (ignoring strand orientation)
for (motif in unique(chm13_gr$motif_id)) {
  mcols(chm13_censat_gr)[[motif]] <- countOverlaps(
    chm13_censat_gr,
    chm13_gr[chm13_gr$motif_id == motif],
    ignore.strand = TRUE
  )
}

# Convert the GRanges object to a dataframe
censat_counts_df <- as.data.frame(mcols(chm13_censat_gr))
motifs <- c("MA0139.2", "MA1929.2", "MA1930.2")

# Calculate CTCF density per CenSat type per motif type
# Density is expressed as number of CTCF sites per 1Mb of genomic sequence
censat_density_df <- censat_counts_df |>
  group_by(type) |>
  summarise(
    n_regions = dplyr::n(),
    total_length = sum(length),
    across(all_of(motifs), sum)
  )

per_chrom_perm_res_files <- fs::dir_ls(
  fs::path(res_path, "misc", "gat", "per_chrom_results"),
  glob = "*.tsv"
)

per_chrom_perm_res_dfs <- purrr::map(
  per_chrom_perm_res_files,
  ~ read_tsv(
    .x,
    show_col_types = FALSE
  )
)

censat_per_chrom_perm_res <- per_chrom_perm_res_dfs |>
  dplyr::bind_rows(.id = "file") |>
  dplyr::mutate(
    seqnames = sub(
      "^chr",
      "",
      stringr::str_extract(file, "chr[0-9XYM]+(?=\\.tsv$)")
    ),
    .keep = "unused"
  ) |>
  dplyr::rename(
    motif = track,
    type = annotation
  ) |>
  dplyr::mutate(sig = case_when(
    qvalue < 0.001 & fold > 1 ~ "***",
    qvalue < 0.01  & fold > 1 ~ "**",
    qvalue < 0.05  & fold > 1 ~ "*",
    .default = " "
  )) |>
  dplyr::select(seqnames, motif, type, qvalue, sig)


gsat_per_chrom_df <- as.data.frame(chm13_censat_gr) |>
    dplyr::filter(type == "gsat")

gsat_density <- gsat_per_chrom_df |>
  dplyr::group_by(seqnames) |>
  dplyr::summarise(
    n_arrays = n(),
    total_len_bp = sum(width),
    across(all_of(motifs), sum),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(all_of(motifs), names_to = "motif", values_to = "sites") |>
  left_join(
      censat_per_chrom_perm_res |>
          dplyr::filter(type == "gsat") |>
          dplyr::mutate(seqnames = paste0("chr", seqnames)),
      by = c("seqnames", "motif")
  ) |>
  dplyr::mutate(
    sig      = tidyr::replace_na(sig, ""),
    density  = sites / total_len_bp * 1e6,
    seqnames = factor(seqnames, levels = paste0("chr", c(1:22, "X", "Y")))
  ) |>
  droplevels()

per_mat_per_chrom_plot <- ggplot(gsat_density, aes(x = seqnames, y = density, fill = motif)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 4.5) +
  facet_wrap(~ motif, ncol = 1, scales = "free_y") +
  labs(x = "Chromosome", y = "CTCF density (sites per Mb)") +
  theme_classic() +
  scale_fill_manual(values = motif_palette, guide = "none") +
  scale_x_discrete(labels = function(x) sub("^chr", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.background = NULL
  )

ggsave(
  fs::path(plots, "per_chrom_per_mat_gsat.pdf"),
  per_mat_per_chrom_plot,
  width = 4,
  height = 3,
  bg = "transparent"
)