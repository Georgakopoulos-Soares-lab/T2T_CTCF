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
    library("ComplexHeatmap", quietly = TRUE, verbose = FALSE)
    library("circlize", quietly = TRUE, verbose = FALSE)
    library("svglite", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata_path for T2T primates
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
    dplyr::rename(name = motif_id) |>
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

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load the IDs for genes present only in the T2T-CHM13 assembly
chm13_only_genes <- read_csv(
  fs::path(data_path, "chm13", "chm13_new_genes.csv"),
  show_col_types = FALSE
)

# Load the GENCODE GFF3 file and convert to dataframe
gencode_gff3 <- import(
  fs::path(data_path, "chm13", "chm13.draft_v2.0.gene_annotation.gff3")
)

gencode_df <- as.data.frame(gencode_gff3)

# Extract gene IDs and filter the gencode dataframe to keep only rows for the known genes
new_gene_ids <- chm13_only_genes$gene_ID
known_gene_rows <- gencode_df |>
    dplyr::filter(type == "gene", !(ID %in% new_gene_ids))

# Convert to GRanges object
known_gene_gr <- makeGRangesFromDataFrame(
    known_gene_rows,
    keep.extra.columns = TRUE
)

# Add promoter annotations
promoter_gr <- promoters(
  known_gene_gr,
  upstream = 2000,
  downstream = 150
)

promoter_df <- as.data.frame(promoter_gr) |>
  mutate(type = "promoter")

# Keep transcripts whose Parent is one of the known genes
transcript_rows <- gencode_df |>
  filter(type %in% c("mRNA", "transcript"), !(Parent %in% new_gene_ids))

transcript_ids <- transcript_rows$ID

# Keep exons/introns/5' and 3' prime UTRs whose
# Parent is one of the known genes
subfeature_rows <- gencode_df |> filter(
  type %in% c("exon", "intron", "five_prime_UTR", "three_prime_UTR"),
  Parent %in% transcript_ids
)

# Combine all and remove strand information
known_gene_features <- bind_rows(subfeature_rows, promoter_df)
known_gene_features$strand <- "*"

# Count how many new genes overlap with CTCF sites
overlaps <- countOverlaps(
  known_gene_gr,
  chm13_gr,
  ignore.strand = TRUE
)

with_ctcf <- sum(overlaps > 0)
without_ctcf <- sum(overlaps == 0)

# Extract the names of known genes that carry CTCF sites
genelist <- as.data.frame(known_gene_gr) |>
  filter(overlaps > 0) |>
  pull(ID)

# Attach overlap counts to the genes so the filter can't get misaligned
known_gene_df <- as.data.frame(known_gene_gr) |>
  mutate(n_ctcf = overlaps)

# Define a named list to rename the different transcript categories for plotting
biotype_labels <- c(
  "protein_coding" = "Protein coding",
  "lncRNA" = "lncRNA",
  "miRNA" = "miRNA",
  "rRNA" = "rRNA",
  "processed_pseudogene" = "Processed pseudogene",
  "unprocessed_pseudogene" = "Unprocessed pseudogene",
  "transcribed_unprocessed_pseudogene" = "Transcribed unprocessed pseudogene",
  "transcribed_processed_pseudogene" = "Transcribed processed pseudogene",
  "StringTie" = "StringTie",
  "TEC" = "TEC",
  "IG_C_gene" = "IG C gene",
  "snRNA" = "snRNA",
  "IG_V_pseudogene" = "IG V pseudogene",
  "rRNA_pseudogene" = "rRNA pseudogene",
  "polymorphic_pseudogene" = "Polymorphic pseudogene",
  "snoRNA" = "snoRNA",
  "misc_RNA" = "Miscellaneous RNA",
  "Mt_tRNA" = "Mitochondrial tRNA",
  "Mt_rRNA" = "Mitochondrial rRNA",
  "IG_V_gene" = "IG V gene",
  "scaRNA" = "scaRNA",
  "TR_J_gene" = "TR J gene",
  "TR_C_gene" = "TR C gene",
  "TR_V_gene" = "TR V gene",
  "unitary_pseudogene" = "Unitary pseudogene",
  "scRNA" = "scRNA",
  "IG_C_pseudogene" = "IG C  pseudogene",
  "unknown_likely_coding" = "Unknown likely  coding",
  "transcribed_unitary_pseudogene" = "Transcribed unitary pseudogene"
)

# Out of the known genes containing at least CTCF we plot the distribution across categories 
known_genes_bar <- known_gene_df |>
  filter(n_ctcf > 0) |>
  count(biotype = gene_biotype) |>
  mutate(
    pct = n / sum(n) * 100,
    biotype = forcats::fct_reorder(biotype, n, .desc = FALSE)
  ) |>
  ggplot(aes(y = biotype, x = pct)) +
  geom_col(width = 0.7, colour = "black", fill = "#636566", linewidth = 0.25) +
  geom_text(aes(label = n), hjust = -0.25, size = 7 / .pt) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_y_discrete(labels = biotype_labels) +
  labs(x = "Percentage (%)", y = NULL) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 6),
    axis.text.y  = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    plot.margin  = margin(5.5, 5.5, 5.5, 5.5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(
  fs::path(plots_path, "known_genes_bar.pdf"),
  known_genes_bar,
  width = 4.13,
  height = 3.10,
  bg = "transparent"
)

# Create a GRanges object from the features of the new genes
known_gene_features_gr <- makeGRangesFromDataFrame(
  known_gene_features,
  keep.extra.columns = TRUE
)

# Find which overlap with CTCF and convert to a dataframe
hits <- findOverlaps(
  chm13_gr,
  known_gene_features_gr,
  ignore.strand = TRUE
)

overlap_df <- data.frame(
  motif = chm13_gr$name[queryHits(hits)],
  feature_type = known_gene_features_gr$type[subjectHits(hits)]
)

# Create a summary of site counts for ploting
overlap_summary <- overlap_df |>
  group_by(motif, feature_type) |>
  summarise(count = n(), .groups = "drop")

# Fix labels for plotting
overlap_summary <- overlap_summary |>
  mutate(
    feature_type = str_to_title(feature_type),
    feature_type = str_replace(feature_type, "Five_prime_utr", "5' UTR"),
    feature_type = str_replace(feature_type, "Three_prime_utr", "3' UTR"),
  )

# CTCF motif palette used throught the project
motif_palette <- c(
  "MA0139.2" = "#A3ABC0",
  "MA1929.2" = "#424D6A",
  "MA1930.2" = "#0E1933"
)

# Create a grouped bar chart with CTCF sites per motif per feature type
# Feature types analyzed: Promoters, exons, introns, 5' UTRs, 3' UTRs
motif_grouped <- ggplot(
  overlap_summary,
  aes(x = count, y = feature_type, fill = motif)
) +
  geom_col(position = "dodge", color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = count, group = motif),
    position = position_dodge(width = 0.9),
    hjust = -0.2,
    size = 6 / .pt
  ) +
  labs(x = "CTCF site count", y = NULL, fill = "motif") +
  scale_x_log10(expand = expansion(mult = c(0, 0.12))) +
  theme_classic() +
  scale_y_discrete(limits = c(
    "3' UTR", "5' UTR", "Intron",
    "Exon", "Promoter"
  ))+
  scale_fill_manual(values = motif_palette) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    plot.background = NULL,
    axis.title.x = element_text(size = 8),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(
  fs::path(plots_path, "known_genes_grouped.pdf"),
  motif_grouped,
  width = 4.13,
  height = 3.10,
  bg = "transparent"
)