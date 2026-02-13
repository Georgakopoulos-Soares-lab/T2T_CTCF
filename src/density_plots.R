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
        library("tidyr", quietly = TRUE, verbose = FALSE)
        library("ggplot2", quietly = TRUE, verbose = FALSE)
        library("ggsci", quietly = TRUE, verbose = FALSE)
        library("patchwork", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

fimo_res = fs::path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data = fs::path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
plots = fs::path(config$plots) # Path to store any generated plots

load_fimo_results <- function(genome, fimo_res_dir) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(fs::path(fimo_res_dir, genome), glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
  ) |>
    select(-c(motif_alt_id, score, `p-value`, `q-value`, matched_sequence)) # Remove unnecessary columns
  )

  res <- dplyr::bind_rows(ctcf_dfs) # Create one complete dataframe with the results

  # Convert the dataframe into a GRanges object for downstream interval manipulation
  ctcf_gr <- GenomicRanges::makeGRangesFromDataFrame(
    res,
    seqnames.field = "sequence_name",
    start.field = "start",
    end.field = "stop",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  
  # Also, create a GRangesList that holds each motif into a separate object
  ctcf_grl <- split(ctcf_gr, ctcf_gr$motif_id)
  
  list(gr = ctcf_gr, grl = ctcf_grl) # Store the two objects (GRanges and GRangesList) 
}

chm13 <- load_fimo_results("chm13", fimo_res)

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

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"), format = "BED")

# Extract CenSat type from the names 
chm13_censat_gr$type <- sub("_.*", "", chm13_censat_gr$name)

# If an HOR has L at the end of the name then this indicates that this is the active alpha-satellite HOR
chm13_censat_gr$type[chm13_censat_gr$type == "hor" & grepl("L\\)$", chm13_censat_gr$name)] <- "hor_l"

# Drop unnecessary collumns
mcols(chm13_censat_gr)[, c("score", "itemRgb", "thick")] <- NULL

# Add a column specifing the region's length
chm13_censat_gr$length <- width(chm13_censat_gr)

# Loop over the CTCF motifs and find how many overlap with the CenSat annotations (ignoring strand orientation)
for (motif in names(chm13$grl)) {
  mcols(chm13_censat_gr)[[motif]] <- countOverlaps(
    chm13_censat_gr,
    chm13$grl[[motif]],
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
    across(motifs, sum)
  )

# Read the chromosome sizes for the CHM13 genome
chm13_chrom_sizes <- read_tsv(
    fs::path(data, "chm13", "chm13v2.0.fa.chrom.sizes"), 
    show_col_types = FALSE,
    col_names = c("chrom", "size"))

# Count genome length by summing the chromosome sizes (including MT)
chm13_genome_size <- sum(chm13_chrom_sizes$size)

# Create an empty dataframe to store the genome-wide densities for each motif
chm13_genome_densities <- data.frame(
  motif = character(),
  density_per_Mb = numeric()
)

# Loop over motifs and calculate the both the CenSat densities and the genome wide densities
for (motif in motifs) {
  censat_density_df[[paste0(motif, "_per_Mb")]] <- censat_density_df[[motif]] / (censat_density_df$total_length / 1e6)
  chm13_genome_densities <- rbind(chm13_genome_densities, data.frame(
    motif = motif,
    density_per_Mb = length(chm13$grl[[motif]]) / (chm13_genome_size / 1e6)
  ))
}

# Define a color palette to be used throughout the plots
motif_palette <- c("MA0139.2" = "#A3ABC0", "MA1929.2" = "#424D6A", "MA1930.2" = "#0E1933")
per_mb_cols <- paste0(motifs, "_per_Mb")

# Pivot the CenSat densities dataframe to be used for plotting
censat_plot_df <- censat_density_df |>
  pivot_longer(
    cols = all_of(per_mb_cols),
    names_to = "motif", 
    values_to = "density_per_Mb") |>
  mutate(motif = gsub("_per_Mb", "", motif))

# Rename the CenSat types for plotting
censat_plot_df$type <- censat_mapping[censat_plot_df$type]

# Create bar plot showing motif densities with vertical lines showing genome-wide density
density_bars <- ggplot(censat_plot_df, aes(x = density_per_Mb, y = type, fill = motif)) +
  geom_vline(
    data = chm13_genome_densities, 
    aes(xintercept = density_per_Mb, color = motif), 
    linetype = "dashed", linewidth = 0.4, alpha = 0.9
  ) +
  geom_col(position = "dodge", color = "black", linewidth = 0.25) +
  labs(x = "Motif density\n(sites per 1Mb)", y = NULL, fill = "Motif") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
        panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
        plot.background = NULL,
        legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8)) +
  scale_fill_manual(values = motif_palette) +
  scale_color_manual(values = motif_palette) +
  scale_y_discrete(limits = c(
    "rDNA array", "Other satellites", "Centromeric transitions",
    "HSat3", "HSat2", "HSat1B", "HSat1A", "γSat", "βSat",
    "αSat divergent", "αSat monomeric", "αSat inactive HOR",
    "αSat active HOR"
  )) +
  guides(fill = guide_legend(), color = "none")

# Convert the previous CenSAt GRanges object containing the number of
# strand-agnostic motif overlaps into a dataframe
chm13_censat_ctcf_overlap_counts <- as.data.frame(chm13_censat_gr)

# Summarize the counts based on CenSat type and chromosome
chm13_censat_ctcf_overlaps_summary <- chm13_censat_ctcf_overlap_counts |>
  group_by(seqnames, type) |>
  summarise(
    # total_MA0139.2 = sum(MA0139.2),
    # total_MA1929.2 = sum(MA1929.2),
    # total_MA1930.2 = sum(MA1930.2),
    total_motifs = sum(MA0139.2 + MA1929.2 + MA1930.2),
    total_length = sum(length),
    .groups = "drop"
  )
# Add density column, which shows the number of overlaps based on total
# region length (summarized across all type instances in a chromosome)
chm13_censat_ctcf_overlaps_summary$density <- chm13_censat_ctcf_overlaps_summary$total_motifs / (chm13_censat_ctcf_overlaps_summary$total_length / 1e6)

# Convert zeros to NA for correct log-scale plotting
chm13_censat_ctcf_overlaps_summary <- chm13_censat_ctcf_overlaps_summary |>
  mutate(density = ifelse(density == 0, NA, density)) |>
  complete(seqnames, type) # Add any missing combinations of type x chromosome
  
# Rename types for plotiing
chm13_censat_ctcf_overlaps_summary$type <- censat_mapping[chm13_censat_ctcf_overlaps_summary$type]

# Remove chr prefix for the plot
chm13_censat_ctcf_overlaps_summary$seqnames <- gsub("^chr", "", chm13_censat_ctcf_overlaps_summary$seqnames)

# Define chromosome order for the plot to be alphanumeric
chr_order <- c(as.character(1:22), "X", "Y")
chm13_censat_ctcf_overlaps_summary$seqnames <- factor(
  chm13_censat_ctcf_overlaps_summary$seqnames,
  levels = chr_order
)

# Heatmap plot showing combined density (all motif types) across each chromosome
density_heatmap <- ggplot(
  chm13_censat_ctcf_overlaps_summary, 
  aes(x = seqnames, y = type, fill = density)) +
  geom_tile(color = "black", linewidth = 0.15) +
  scale_fill_viridis_c(
    option="G", 
    trans = "log10", 
    na.value = "grey95", 
    direction = -1,
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    plot.background = NULL,
    axis.title.x = element_text(size = 8),
    legend.box.just = "center",
    legend.justification = "center",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.4, "cm"),
    legend.key.height = unit(0.4, "cm")
  ) +
  labs(x = "Chromosome", y = NULL, fill = "CTCF density\n     (log10)") +
  scale_y_discrete(limits = c(
    "rDNA array", "Other satellites", "Centromeric transitions",
    "HSat3", "HSat2", "HSat1B", "HSat1A", "γSat", "βSat",
    "αSat divergent", "αSat monomeric", "αSat inactive HOR",
    "αSat active HOR"
  ))

# Using patchwork combine the two density plots into one and save
combined_density_plots <- (density_bars | density_heatmap) + 
  plot_layout(widths = c(1, 3)) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.2, "cm")
  )

ggsave(
  fs::path(plots, "combined_densities.svg"),
  combined_density_plots, 
  width = 5.34,
  height = 3.4
)