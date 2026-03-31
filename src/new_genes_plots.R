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

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results
fimo_res <- fs::path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data <- fs::path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
plots <- fs::path(config$plots) # Path to store any generated plots

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

chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load the IDs for genes present only in the T2T-CHM13 assembly
chm13_only_genes <- read_csv(
    fs::path(data, "chm13", "chm13_new_genes.csv"),
    show_col_types = FALSE
)

# Load the GENCODE GFF3 file and convert to dataframe
gencode_gff3 <- import(fs::path(data, "chm13", "chm13.draft_v2.0.gene_annotation.gff3"))
gencode_df <- as.data.frame(gencode_gff3)

# Extract gene IDs and the corresponding rows from the GENCODE file
gene_ids <- chm13_only_genes$gene_ID
gene_rows <- gencode_df |> filter(type == "gene", ID %in% gene_ids)

# Convert to GRanges object
gene_gr <- makeGRangesFromDataFrame(gene_rows, keep.extra.columns = TRUE)

# Add promoter annotations
promoter_gr <- promoters(gene_gr, upstream = 2000, downstream = 150)

promoter_df <- as.data.frame(promoter_gr) |>
  mutate(type = "promoter")

# Keep transcripts whose Parent is one of the new T2T-CHM13-only genes
transcript_rows <- gencode_df |> filter(type %in% c("mRNA", "transcript"), Parent %in% gene_ids)
transcript_ids <- transcript_rows$ID

# Keep exons/introns/5' and 3' prime UTRs whose Parent is one of the new T2T-CHM13-only genes
subfeature_rows <- gencode_df |> filter(
    type %in% c("exon", "intron", "five_prime_UTR", "three_prime_UTR"),
    Parent %in% transcript_ids
)

# Combine all and remove strand information
new_gene_features <- bind_rows(subfeature_rows, promoter_df)
new_gene_features$strand <- "*"

# Count how many new genes overlap with CTCF sites
overlaps <- countOverlaps(gene_gr, chm13_gr, ignore.strand = TRUE)
with_ctcf <- sum(overlaps > 0)
without_ctcf <- sum(overlaps == 0)

# Extract the names of T2T-CHM13 novel genes that carry CTCF sites
genelist <- as.data.frame(gene_gr) |> 
  filter(overlaps > 0) |>
  pull(ID)

# Create a stacked bar chart showing the distriburion of gene categories
# for those genes
new_genes_stacked <- chm13_only_genes |>
  filter(gene_ID %in% genelist) |>
  count(biotype) |>
  mutate(pct = n / sum(n) * 100) |>
  ggplot(aes(y = "", x = pct, fill = biotype)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.25) +
  labs(
    y = NULL,
    x = "Percentage (%)",
    fill = "Biotype",
  ) +
  theme_classic() +
  scale_fill_manual(values = c(
    "protein_coding" = "#2B5C8A",
    "lncRNA" = "#8B7BAA",
    "miRNA" = "#E8A838",
    "rRNA" = "#D4884A",
    "processed_pseudogene" = "#1B3A5C",
    "unprocessed_pseudogene" = "#3D6E99",
    "transcribed_unprocessed_pseudogene" = "#6A9BC5",
    "transcribed_processed_pseudogene" = "#A3C4DE",
    "StringTie" = "#5B5B5B",
    "TEC" = "#C2A8D1",
    "IG_C_gene" = "#4A7A6F"
  )) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 6),
  ) +
  scale_y_discrete(expand = expansion(add = 0.1))

ggsave(
    fs::path(plots, "new_genes_stacked.pdf"),
    new_genes_stacked,
    width = 3,
    height = 0.75
)

# Create a dataframe with the values of genes with CTCF vs those with no CTCF sites
pie_df <- data.frame(
  category = c("With CTCF", "Without CTCF"),
  count    = c(with_ctcf, without_ctcf)
)

# Calculate percentages
pie_df$percentage <- round(pie_df$count / sum(pie_df$count) * 100, 2)
pie_df$label <- paste0(pie_df$category, "\n", pie_df$count, " (", pie_df$percentage, "%)")

# Plot pie chart and save 
pie <- ggplot(pie_df, aes(x = "", y = count, fill = category)) +
  geom_col(width = 10, color= "black") +
  coord_polar(theta = "y", start = 2 * pi/3 + 0.2) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2
  ) +
  scale_fill_manual(values = c("With CTCF" = "#A3ABC0", "Without CTCF" = "#CD979A")) +
  ggtitle("Novel T2T-CHM13 genes\nwith CTCF sites") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold")
  )

ggsave(
    fs::path(plots, "new_genes_pie.pdf"),
    pie,
    width = 2,
    height = 2
)

# Create a GRanges object from the features of the new genes
new_gene_features_gr <- makeGRangesFromDataFrame(new_gene_features, keep.extra.columns = TRUE)

# Find which overlap with CTCF and convert to a dataframe
hits <- findOverlaps(chm13_gr, new_gene_features_gr, ignore.strand = TRUE)
overlap_df <- data.frame(
  motif    = chm13_gr$motif_id[queryHits(hits)],
  feature_type  = new_gene_features_gr$type[subjectHits(hits)]
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
motif_palette <- c("MA0139.2" = "#A3ABC0", "MA1929.2" = "#424D6A", "MA1930.2" = "#0E1933")

# Create a grouped bar chart with CTCF sites per motif per feature type
# Feature types analyzed: Promoters, exons, introns, 5' UTRs, 3' UTRs
motif_grouped <- ggplot(overlap_summary, aes(x = count, y = feature_type, fill = motif)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.25) +
  labs(x = "CTCF Site Count", y = NULL, fill = "motif") +
  scale_x_log10() +
  theme_classic()+ 
  scale_y_discrete(limits = c(
    "3' UTR", "5' UTR", "Intron",
    "Exon", "Promoter"
  ))+
  scale_fill_manual(values = motif_palette) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.position = "none",
    plot.background = NULL,
    axis.title.x = element_text(size = 8),
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
  )

ggsave(
    fs::path(plots, "new_genes_grouped.pdf"),
    motif_grouped,
    width = 3,
    height = 2.25
)