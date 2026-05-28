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
    library("Biostrings", quietly = TRUE, verbose = FALSE)
    library("Rsamtools", quietly = TRUE, verbose = FALSE)
    library("plyranges", quietly = TRUE, verbose = FALSE)
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

# Load CHM13 chromosome sizes
chm13_chrom_sizes <- read_tsv(
  fs::path(
    data, "chm13", "chm13v2.0.chrom.sizes"
  ),
  show_col_types = FALSE,
  col_names = c("chrom", "size")
)

# Calculate CTCF genome-wide density across all 3 motifs
genome_density <- (length(chm13_gr) / sum(chm13_chrom_sizes$size)) * 1e6

# Load CHM13 pseudohomologous regions (PHRs)
# present in the acrocentric short arms
chm13_phrs <- rtracklayer::import(
  fs::path(
    data, "chm13", "chm13_acro_PHRs.bed"
  )
)

# Find how many CTCF sites overlap with each PHR
chm13_phrs$ctcf_sites <- countOverlaps(
  chm13_phrs,
  chm13_gr,
  ignore.strand = TRUE
)

# Define the path where the T2T-CHM13v2.0 genome FASTA is store
fasta_file <- fs::path(data, "chm13", "chm13v2.0.fa.gz")

# Provide the FASTA file as a pointer for R to use
# Make sure the genome file is indexed using Samtools
fa <- Rsamtools::FaFile(fasta_file)

# Get the underlying sequence for each window
chm13_phr_sequences <- Biostrings::getSeq(fa, chm13_phrs)

# Calculate GC content (counts / region's length)
# and add it as metadata to the windows GRanges
gc_count <- as.integer(
  Biostrings::letterFrequency(
    chm13_phr_sequences, "GC"
  )
)

chm13_phrs <- chm13_phrs |>
  mutate(
    gc_count = gc_count,
    gc_content = gc_count / width(chm13_phrs)
  )

# Convert to dataframe
chm13_phrs_df <- as.data.frame(chm13_phrs)

# In the data the name column is the Shannon Diversity Index (SDI)
# SDI shows the positional homology entropy, and it serves to show
# regions where contigs can be described as mosaics of different
# reference chromosomes
chm13_phrs_df$sdi <- as.numeric(chm13_phrs_df$name)

# Calculate density of CTCF sites in each PHR
chm13_phrs_df <- chm13_phrs_df |>
  dplyr::select(-c(strand, name)) |>
  dplyr::mutate(
    density = (ctcf_sites / width) * 1e6,
    width_mb = width / 1e6
  )

# Spearman correlation between CTCF density and SDI
sp <- cor.test(
  chm13_phrs_df$density,
  chm13_phrs_df$sdi,
  method = "spearman",
  exact = FALSE
)

# Format Spearman p-value
p_label <- if (sp$p.value < 1e-3) {
  "italic(p) < 0.001"
} else {
  sprintf("italic(p) == %.3f", sp$p.value)
}

# Create label for plotting
cor_label <- sprintf(
  '"Spearman"~rho == %.2f*","~%s',
  sp$estimate,
  p_label
)

# Create a plot showing the relationship between SDI and CTCF density
sdi_corr <- ggplot(
  chm13_phrs_df,
  aes(x = sdi, y = density)
) +
  geom_hline(
    yintercept = genome_density,
    linetype = "dashed",
    linewidth = 0.5,
    color = "#383d50"
  ) +
  # Fit a LOESS curve to the data to show the relationship
  geom_smooth(
    method = "loess",
    span = 1,
    se = TRUE,
    color = NA,
    fill = "grey75"
  ) +
  geom_point(
    data = ~ dplyr::arrange(.x, dplyr::desc(width_mb)),
    aes(size = width_mb, fill = score),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  # Fit a LOESS curve to the data to show the relationship
  geom_smooth(
    method = "loess",
    span = 1,
    se = TRUE,
    color = "grey25",
    fill = NA
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_label,
    parse = TRUE,
    hjust = -0.05,
    vjust = 1.4,
    size = 2.2
  ) +
  scale_size_continuous(
    range = c(1, 5)
  ) +
  scale_fill_viridis_c(
    option = "cividis",
    alpha = 0.7
  ) +
  coord_cartesian(ylim = c(0, 115)) +
  labs(
    x = "Shannon diversity Index",
    y = "CTCF density (sites / 1Mb)",
    fill = "Supporting\nhaplotypes",
    size = "Region size (Mb)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "none",
    # legend.position = "bottom",
    # legend.justification = "center",
    # legend.title = element_text(size = 6),
    # legend.text = element_text(size = 6),
    # legend.key.width = unit(0.45, "cm"),
    # legend.key.height = unit(0.2, "cm")
  )

ggsave(
  fs::path(plots, "acro_ctcf_sdi_correlation.pdf"),
  sdi_corr,
  width = 2.1,
  height = 2,
  bg = "transparent"
)