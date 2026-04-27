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
    library("BiocParallel", quietly = TRUE, verbose = FALSE)
    library("ggridges", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

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

  # Make the results strand-agnostic
  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# A mapping that goes from abbreviations to full names (used for plotting)
censat_mapping <- c(
  "ct" = "Centromeric\ntransitions",
  "censat" = "Other satellites",
  "gsat" = "γSat",
  "bsat" = "βSat",
  "hsat3" = "HSat3",
  "hsat2" = "HSat2",
  "hsat1b" = "HSat1B",
  "hsat1a" = "HSat1A",
  "mon" = "αSat monomeric",
  "dhor" = "αSat divergent",
  "hor" = "αSat inactive HOR",
  "hor_l" = "αSat active HOR",
  "rdna" = "rDNA array"
)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# If an HOR has L at the end of the name then
# this indicates that this is the active alpha-satellite HOR
chm13_censat_gr$type[chm13_censat_gr$type == "hor" & grepl("L\\)$", chm13_censat_gr$name)] <- "hor_l"

# Drop unnecessary collumns
mcols(chm13_censat_gr)[, c("score", "itemRgb", "thick")] <- NULL

# Store the BigWig file containg conservation
# scores among primates as a pointer for fast access
conserv_bw <- BigWigFile(fs::path(data, "primates", "approach2.scores_v0.3.bw"))

# Create a function which adds mean conservation scores
# to a GRanges object
add_mean_conservation <- function(
  gr,
  conserv_bw
) {

  # Sort to improve BigWig I/O locality
  gr <- sort(gr)

  # Preserve original order so we can map scores back correctly
  gr$.orig_idx <- seq_along(gr)

  # Split by chromosome for parallel processing
  by_chr <- split(gr, GenomeInfoDb::seqnames(gr))
  by_chr <- by_chr[lengths(by_chr) > 0]

  scores_list <- BiocParallel::bplapply(
    by_chr,
    function(chr_gr) {
      s <- summary(conserv_bw, which = chr_gr, type = "mean")
      idx <- chr_gr$.orig_idx
      vals <- vapply(s, function(x) x$score[1], numeric(1))
      data.frame(idx = idx, score = vals)
    },
    BPPARAM = BiocParallel::SnowParam(workers = 6)
  )

  combined <- do.call(rbind, scores_list)
  scores <- combined$score[order(combined$idx)]

  gr <- gr[order(gr$.orig_idx)]
  gr$.orig_idx <- NULL
  gr$mean_conservation <- scores

  gr
}

# Add conservation scores to both FIMO results and
# CenSat regions
chm13_censat_gr <- add_mean_conservation(
  chm13_censat_gr,
  conserv_bw
)

chm13_gr <- add_mean_conservation(
  chm13_gr,
  conserv_bw
)

# Find overlaps between CTCF sites and CenSat regions
hits <- findOverlaps(chm13_gr, chm13_censat_gr, ignore.strand = TRUE)

# Initialize the censat_type column as NA for all CTCF sites
chm13_gr$censat_type <- NA

# For each CTCF site that overlaps a censat region, assign the type
overlap_types <- split(
  chm13_censat_gr$type[subjectHits(hits)],
  queryHits(hits)
)

# Assign back to the CTCF sites and
# map the short codes to full CenSat names
collapsed <- vapply(overlap_types, paste, character(1), collapse = ";")

chm13_gr$censat_type <- NA_character_
chm13_gr$censat_type[as.integer(names(collapsed))] <- collapsed

# Convert the GRanges object to a dataframe for easier plotting
ctcf_sites_df <- as.data.frame(chm13_gr)

# Add the full names to each CenSat type and drop sites that do not overlap
# with CenSat types
ctcf_sites_df$censat_type <- censat_mapping[ctcf_sites_df$censat_type]
ctcf_sites_df <- ctcf_sites_df |> filter(censat_type != "NA")

# Count how many CTCF sites exist per CenSat type
ctcf_counts <- ctcf_sites_df |>
  dplyr::mutate(
    censat_type = factor(censat_type, levels = unname(censat_mapping))
  ) |>
  dplyr::count(censat_type, .drop = FALSE) |>
  tidyr::complete(censat_type = unname(censat_mapping), fill = list(n = 0)) |>
  dplyr::mutate(label = paste0("n=", n))

# Add full names to each CenSat type
chm13_censat_gr$censat_type <- censat_mapping[chm13_censat_gr$type]

# Combine conservation data for the full CenSat regions
# and the CTCF sites within them
plot_df <- rbind(
  ctcf_sites_df |>
    dplyr::select(
      c("censat_type", "mean_conservation")
    ) |>
    dplyr::mutate(cat = "CTCF"),
  as.data.frame(chm13_censat_gr) |>
    dplyr::select(
      c("censat_type", "mean_conservation")
    ) |>
    dplyr::mutate(cat = "CenSat")
)

# Order the categories so that the CTCF density line
# overlays the CenSat background
plot_df$cat <- factor(plot_df$cat, levels = c("CenSat", "CTCF"))

# This is obtained from UCSC's bigWigInfo
genome_mean_conservation <- 0.224605

# Perform Kolmogorov-Smirnoff test to check if the distributions of
# conservation scores differ between CenSat regions and CTCF sites
ks_results <- plot_df |>
  group_by(censat_type) |>
  filter(any(cat == "CTCF") & any(cat == "CenSat")) |>
  summarise(
    test = list(ks.test(
      mean_conservation[cat == "CTCF"],
      mean_conservation[cat == "CenSat"]
    )),
    n_ctcf = sum(cat == "CTCF"),
    n_censat = sum(cat == "CenSat"),
    .groups = "drop"
  ) |>
  mutate(
    D = map_dbl(test, "statistic"),
    p_value = map_dbl(test, "p.value"),
    p_adj = p.adjust(p_value, method = "BH")
  ) |>
  select(-test) |>
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Create a ridgeline plot showing the density of points based on conservation scores
censat_conserv <- ggplot(
  plot_df,
  aes(
    x = mean_conservation,
    y = censat_type,
    fill = cat
  )
) +
  geom_vline(
    xintercept = genome_mean_conservation,
    linetype = "dashed",
    linewidth = 0.75,
    color = "#61573D"
  ) +
  geom_vline(
    xintercept = 0.75,
    linetype = "dashed",
    linewidth = 0.75,
    color = "#e1cc97"
  ) +
  geom_density_ridges(
    color = "black",
    alpha = 0.6,
    scale = 1,
    from = 0,
    to = 1,
    linewidth = 0.75,
    rel_min_height = 0.01
  ) +
  geom_text(
    data = ctcf_counts,
    aes(x = 0.75, y = censat_type, label = label),
    hjust = 0, vjust = -0.6,
    size = 2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ks_results,
    aes(x = 0.75, y = censat_type, label = sig),
    hjust = 0, vjust = -1,
    size = 2.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Mean conservation",
    y = NULL,
    title = NULL
  ) +
  scale_y_discrete(limits = c(
    "rDNA array", "Other satellites", "Centromeric\ntransitions",
    "HSat3", "HSat2", "HSat1B", "HSat1A", "γSat", "βSat",
    "αSat divergent", "αSat monomeric", "αSat inactive HOR",
    "αSat active HOR"
  ),
  expand = expansion(add = c(0.1, 1))
  ) +
  scale_fill_manual(
    values = c(
      "CTCF" = "#616883",
      "CenSat" = "gray70"
    )
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    plot.background = NULL
  )

ggsave(
  fs::path(plots, "primate_censat_conserv.svg"),
  censat_conserv,
  width = 2.5,
  height = 4,
  bg = "transparent"
)