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
    library("patchwork", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
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

# Keep a list with the three CTCF motifs we analyze
motifs <- unique(chm13_gr$motif_id)

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

gsat_counts_df <- censat_counts_df |>
  dplyr::filter(type == "gsat") # Keep only gamma-satellites

gsat_counts_df$family <- str_extract(
  gsat_counts_df$name, "(?<=\\().*?(?=\\))"
) # Extract gamma-satellite family

# A mapping that normalizes the names of the different gamma-satellite types
gsat_mapping <- c(
  "GSATII" = "γSatII",
  "GSATII/TAR" = "γSatII-TAR composite",
  "GSATX;GSATII" = "γSatII-γSatX composite",
  "GSAT" = "γSat",
  "GSATII,TAR" = "γSatII-TAR composite",
  "TAR/GSATII" = "γSatII-TAR composite",
  "GSATX" = "γSatX",
  "GSATII,rnd-6_family-1431" = "γSatII",
  "GSATII/TAR1" = "γSatII-TAR composite"
)

# Rename the gamma-satellite families for plotting
gsat_counts_df$family <- gsat_mapping[gsat_counts_df$family]

# Calculate CTCF density per gamma-satellite family per motif type
# Density is expressed as number of CTCF sites per 1Mb of genomic sequence
gsat_density <- gsat_counts_df |>
  group_by(family) |>
  summarise(
    n_regions = dplyr::n(),
    total_length = sum(length),
    across(all_of(motifs), sum)
  )

# Read the chromosome sizes for the CHM13 genome
chm13_chrom_sizes <- read_tsv(
  fs::path(data, "chm13", "chm13v2.0.chrom.sizes"),
  show_col_types = FALSE,
  col_names = c("chrom", "size")
)

# Count genome length by summing the chromosome sizes (including MT)
chm13_genome_size <- sum(chm13_chrom_sizes$size)

# Create an empty dataframe to store the genome-wide densities for each motif
chm13_genome_densities <- data.frame(
  motif = character(),
  density_per_Mb = numeric()
)

# Loop over motifs and calculate both
# the CenSat densities and the genome wide densities
for (motif in motifs) {
  gsat_density[[paste0(motif, "_per_Mb")]] <- gsat_density[[motif]] / (gsat_density$total_length / 1e6)
  chm13_genome_densities <- rbind(chm13_genome_densities, data.frame(
    motif = motif,
    density_per_Mb = length(chm13_gr[chm13_gr$motif_id == motif]) / (chm13_genome_size / 1e6)
  ))
}

# Define a palette for plotting
motif_palette <- c(
  "MA0139.2" = "#A3ABC0",
  "MA1929.2" = "#424D6A",
  "MA1930.2" = "#0E1933"
)
per_mb_cols <- paste0(motifs, "_per_Mb") # Density columns

# Load the results of the permutation test we performed using GAT
gsat_gat_res <- read_tsv(
  fs::path(res_path, "misc", "gat", "GAT_results_gsat.tsv"),
  show_col_types = FALSE
) |>
  dplyr::rename(motif = track, family = annotation) |>
  dplyr::mutate(sig = case_when(
    qvalue < 0.001 & fold > 1 ~ "***",
    qvalue < 0.01  & fold > 1 ~ "**",
    qvalue < 0.05  & fold > 1 ~ "*",
    .default = ""
  )) |>
  dplyr::select(motif, family, qvalue, sig)

# A mapping that normalizes the names of the different gamma-satellite types
gsat_gat_mapping <- c(
  "GSATX" = "γSatX",
  "GSAT" = "γSat",
  "GSATII" = "γSatII",
  "GSATII-TAR composite" = "γSatII-TAR composite",
  "GSATII-GSATX composite" = "γSatII-γSatX composite"
)

gsat_gat_res$family <- gsat_gat_mapping[gsat_gat_res$family]

# Pivot the gamma-satellite density dataframe to be easier for plotting
gsat_plot_df <- gsat_density |>
  pivot_longer(
    cols = all_of(per_mb_cols),
    names_to = "motif",
    values_to = "density_per_Mb"
  ) |>
  mutate(motif = gsub("_per_Mb", "", motif)
  ) |>
  dplyr::left_join(
    gsat_gat_res, by = c("motif", "family")
  )

# Create a barplot showing the density per gamma-satellite family
gsat_density_barplot <- ggplot(
  gsat_plot_df,
  aes(x = density_per_Mb, y = family, fill = motif)
) +
  geom_vline(
    data = chm13_genome_densities,
    aes(xintercept = density_per_Mb, color = motif),
    linetype = "dashed", linewidth = 0.4, alpha = 0.9
  ) +
  geom_col(
    position = "dodge",
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = sig, group = motif),
    position = position_dodge(width = 0.9),
    vjust = 0.8,
    hjust = -0.2, # push slightly to the right of the bar
    size = 3
  ) +
  labs(x = "Motif density (sites per 1Mb)", y = NULL, fill = "Motif") +
  theme_classic() +
  scale_fill_manual(values = motif_palette) +
  scale_color_manual(values = motif_palette) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_y_discrete(limits = c(
    "γSat", "γSatII", "γSatII-TAR composite", "γSatII-γSatX composite", "γSatX"
  )) +
  guides(fill = guide_legend(), color = "none") +
  theme(panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
        plot.background = NULL,
        legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8))

# Create GRanges object with only the gamma-satellite regions
gsat_gr <- chm13_censat_gr[chm13_censat_gr$type == 'gsat']

# We find the exact CTCF sites that overlap with gamma-satellite regions
# in a strand-agnostic manner
hits <- findOverlaps(chm13_gr, gsat_gr, ignore.strand=TRUE)
ctcf_gsat_overlaps <- chm13_gr[queryHits(hits)]

# We add back the gamma-satellite information for each overlap
mcols(ctcf_gsat_overlaps) <- cbind(mcols(ctcf_gsat_overlaps), mcols(gsat_gr[subjectHits(hits)]))

# We create a pointer to the BigWigFile containing long-read methylation values
# obtained with nanopolish for the CHM13 cell line
methylation <- BigWigFile(
  fs::path(
    data,
    "chm13",
    "chm13v2.0_CHM13_CpG_ont_guppy3.6.0_nanopolish0.13.2.bw"
  )
)

# For each CTCF site that overlaps with gamma-satellite regions we get the mean methylation
# of the interval and we add back this value to the original dataframe
suppressWarnings({
  mean_meth_scores <- summary(
    methylation,
    which = ctcf_gsat_overlaps,
    type = "mean"
  )
  ctcf_gsat_overlaps$mean_meth <- sapply(
    mean_meth_scores,
    function(x) x$score[1]
  )
  ctcf_gsat_overlaps_df <- as.data.frame(ctcf_gsat_overlaps)
  ctcf_gsat_overlaps_df$family <- str_extract(
    ctcf_gsat_overlaps_df$name, "(?<=\\().*?(?=\\))"
  )
  # Rename the gamma-satellite families for plotting
  ctcf_gsat_overlaps_df$family <- gsat_mapping[ctcf_gsat_overlaps_df$family]
})

# We create a boxplot showing the distribution of methylation scores across the
# CTCF sites that overlap with gamma-satellite regions
ctcf_gsat_meth_boxplot <- ctcf_gsat_overlaps_df |>
  ggplot(aes(x = mean_meth, y = family, fill = motif_id)) +
  geom_vline(
    xintercept = 0.3, color = "#247cd3",
    linetype = "dashed", linewidth = 0.75
  ) +
  geom_vline(
    xintercept = 0.7, color = "#ca2e40",
    linetype = "dashed", linewidth = 0.75
  ) +
  geom_boxplot(
    color = "black", linewidth = 0.3,
    outlier.shape = NA, alpha = 0.85
  ) +
  geom_jitter(
    aes(color = motif_id), size = 0.6,
    alpha = 1, color = "black",
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  ) +
  # scale_color_manual(values = motif_palette) +
  labs(x = "Motif mean methylation", y = NULL, fill = "Motif") +
  scale_fill_manual(values = motif_palette) +
  guides(fill = guide_legend(), color = "none") +
  # Order the y-axis ticks the same as the bar plot
  scale_y_discrete(limits = c(
    "γSat", "γSatII", "γSatII-TAR composite", "γSatII-γSatX composite", "γSatX"
  )) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
        plot.background = NULL,
        legend.position = "none",
        # remove them for space since we are creating a combined plot
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8))

# We create a combined figure for density and methylation
combined_density_plots <- (gsat_density_barplot | ctcf_gsat_meth_boxplot) +
  plot_layout(widths = c(1, 1.5)) &
  theme(
    legend.position = "none",
    legend.justification = "center",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.2, "cm"),
    plot.background = NULL,
  )

ggsave(
  fs::path(plots, "gsat_plus_meth.svg"),
  combined_density_plots,
  width = 5.8,
  height = 2.5,
  bg = "transparent"
)