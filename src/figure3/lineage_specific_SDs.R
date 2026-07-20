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
    library("patchwork", quietly = TRUE, verbose = FALSE)
    library("EnrichedHeatmap", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
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

  # Make the object strand-agnostic
  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

lineage_specific_SDs <- readr::read_tsv(
  fs::path(data_path, "primates", "lineage_specific_SDs.txt"),
  show_col_types = FALSE,
  col_names = c(
    "species", "seqnames", "start",
    "end", "sd_status", "lineage", "haplotype_variability",
    "overlapping_novel_transcripts", "rearrangement_breakpoints"
  ),
  skip = 2
)

primate_name_map <- c(
  "GGO" = "Gorilla gorilla gorilla",
  "PAB" = "Pongo abelii",
  "PPA" = "Pan paniscus",
  "PPY" = "Pongo pygmaeus",
  "PTR" = "Pan troglodytes",
  "SSY" = "Symphalangus syndactylus",
  "GorGor" = "Gorilla gorilla gorilla",
  "PonAbe" = "Pongo abelii",
  "PanPan" = "Pan paniscus",
  "PonPyg" = "Pongo pygmaeus",
  "PanTro" = "Pan troglodytes",
  "SymSyn" = "Symphalangus syndactylus"
)

lineage_specific_SDs$species <- primate_name_map[lineage_specific_SDs$species]

lineage_specific_SDs_gr <- makeGRangesFromDataFrame(
  lineage_specific_SDs,
  keep.extra.columns = TRUE
)

primate_paths <- fs::dir_ls(
  fs::path(fimo_res, "primates"),
  type = "directory"
)

names(primate_paths) <- fs::path_file(primate_paths)

gr_list <- lapply(primate_paths, load_fimo_results)

gr_list <- Map(
  function(gr, sp) {
    gr$species <- sp
    gr
  },
  gr_list,
  names(gr_list)
)

suppressWarnings({
  primate_ctcf_gr <- do.call(c, unname(gr_list))
})

load_assembly_report <- function(
  assembly_report_path
) {
  # Load the NCBI assembly report containing the map
  # to go from RefSeq accession to chromosome name
  assembly_report <- readr::read_tsv(
    assembly_report_path,
    comment = "#",
    col_names = c(
      "chromosome", "assigned_role", "assigned_molecule",
      "assigned_molecule_location", "genBank_accn", "relationship",
      "refseq_accn", "assembly_unit", "length", "UCSC_name"
    ),
    show_col_types = FALSE
  )
  assembly_report <- assembly_report |> filter(`assigned_role` == 'assembled-molecule')
  assembly_report
}
# Initialize an empty dataframe instead of a list
primate_genome_densities <- data.frame() 

# Loop over the different directories in the FIMO results path
# We need only the CHM13 and primates folders
for (directory in dir_ls(fimo_res, type = "directory")) {
  if (path_file(directory) == "primates") {
    for (primate_dir in dir_ls(directory, type = "directory")) {
      species <- path_file(primate_dir)

      # Load the NCBI assembly report
      assembly_report <- load_assembly_report(
        fs::dir_ls(
          fs::path(data_path, "primates", species),
          glob = "*_pri_assembly_report.txt"
        )
      )

      fimo_gr <- primate_ctcf_gr[primate_ctcf_gr$species == species]

      # Calculate density and store in a temporary single-row dataframe
      current_species_data <- as.data.frame(fimo_gr) |>
        dplyr::count() |>
        dplyr::mutate(
            species = species, # Add the species name as a column
            density = (n / sum(assembly_report |> 
                                 filter(genBank_accn != "na") |> 
                                 pull(length))) * 1e6
        )
      
      # Append the current species data to the main dataframe
      primate_genome_densities <- dplyr::bind_rows(primate_genome_densities, current_species_data)
    }
  }
}

primate_genome_densities$species <- primate_name_map[primate_genome_densities$species]

primate_ctcf_gr$species <- primate_name_map[primate_ctcf_gr$species]

lineage_specific_SDs_gr$overlaps <- 0L

count_by_species <- function(query, subject) {
  out <- integer(length(query))
  for (sp in unique(query$species)) {
    idx <- which(query$species == sp)
    out[idx] <- countOverlaps(query[idx],
                              subject[subject$species == sp],
                              ignore.strand = TRUE)
  }
  out
}

suppressWarnings({
  lineage_specific_SDs_gr$overlaps <- count_by_species(lineage_specific_SDs_gr, primate_ctcf_gr)
})

plot_df <- as.data.frame(lineage_specific_SDs_gr)

plot_df <- plot_df |>
  dplyr::mutate(
    density = ((overlaps / width) * 1e6)
  ) |>
  dplyr::filter(
    !sd_status %in% c("unmapped", NA),
    sd_status != "SharedByAll"
  )

category_map = c(
  "Expanded" = "Expanded",
  "Expanded&New" = "Expanded\nand new",
  "New" = "New"
)

plot_df$sd_status <- category_map[plot_df$sd_status]

species_order <- c(
  "Symphalangus syndactylus",
  "Pongo abelii", 
  "Pongo pygmaeus",
  "Gorilla gorilla gorilla",
  "Pan paniscus",
  "Pan troglodytes"
)

plot_df$species <- factor(plot_df$species, levels = species_order)
primate_genome_densities$species <- factor(primate_genome_densities$species, levels = species_order)

species_palette <- c(
  "Symphalangus syndactylus" = "#846081",
  "Pongo abelii" = "#B9CFDB",
  "Pongo pygmaeus" = "#D9CAB6",
  "Gorilla gorilla gorilla" = "#587A88",
  "Pan paniscus" = "#C08D66",
  "Pan troglodytes" = "#8D9A73"
)

sd_violin <- ggplot(plot_df, aes(y = species, x = density, fill = species)) +
  geom_violin(
    data = plot_df,
    trim = TRUE,
    scale = "width",
    color = "black",
    linewidth = 0.25,
    alpha = 1
  ) +
  geom_jitter(
    data = plot_df,
    aes(x = density, y = species, fill = species),
    shape = 21,
    height = 0.15,
    width = 0,
    color = "gray15",
    stroke = 0.1,
    size = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = scales::alpha(species_palette, 0.7)
  ) +
  geom_point(
    data = primate_genome_densities,
    aes(x = density, y = species),
    color = "#9d2658",
    shape = 4,
    size = 1.5,
    stroke = 0.6,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ sd_status, nrow = 1) +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 10, 100, 1000, 10000),
    labels = c(0, 10, 100, 1000, 10000)
  ) +
  labs(
    y = NULL, x = "CTCF density") +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray70", linewidth = 0.3),
    plot.background = NULL,
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.position = "none"
  )

ggsave(
  fs::path(plots, "primate_SD_violin.pdf"),
  sd_violin,
  bg = "transparent",
  height = 2,
  width = 5
)