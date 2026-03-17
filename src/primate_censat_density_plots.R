suppressPackageStartupMessages({
  suppressWarnings({
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("ggtree", quietly = TRUE, verbose = FALSE)
    library("ape", quietly = TRUE, verbose = FALSE)
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("GenomeInfoDb", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores the
# paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

# Path containing FIMO results for CHM13, HG002,
# T2T primates, and human pangenome assemblies
fimo_res <- fs::path(config$fimo_res)

# Path containing metadata for CHM13, HG002,
# T2T primates, and human pangenome assemblies
data <- fs::path(config$data)

plots <- fs::path(config$plots) # Path to store any generated plots

# Load primate & human phylogenetic tree
tree <- ape::read.tree(
  fs::path(data, "multiple_sequence_alignments", "8-t2t-apes-2023v2.nh")
)

# Remove hg38 from the plotting
tree <- drop.tip(tree, c("hg38"))

# Create a map to convert GenBank accession names to species names
name_map <- c(
  "GCA_028885655.2" = "Pongo abelii",
  "GCA_028885625.2" = "Pongo pygmaeus",
  "GCA_029281585.2" = "Gorilla gorilla gorilla",
  "GCA_029289425.2" = "Pan paniscus",
  "GCA_028858775.2" = "Pan troglodytes",
  "GCA_028878055.2" = "Symphalangus syndactylus",
  "hs1" = "Homo sapiens"
)

# Rename the GenBank accessions
tree$tip.label <- name_map[tree$tip.label]

load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
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

  # Select only the needed columns
  assembly_report <- assembly_report |> dplyr::select(c(
    "chromosome", "refseq_accn"
  ))

  assembly_report
}

# Create an empty named list to store the FIMO GRanges objects
# across species
fimo_grs <- list()

# Loop over the different directories in the FIMO results path
# We need only the CHM13 and primates folders
for (directory in dir_ls(fimo_res, type = "directory")) {
  if (path_file(directory) == "chm13") {

    species <- path_file(directory)
    fimo_grs[[species]] <- load_fimo_results(directory)

  } else if (path_file(directory) == "primates") {
    for (primate_dir in dir_ls(directory, type = "directory")) {
      species <- path_file(primate_dir)

      # Load the NCBI assembly report
      assembly_report <- load_assembly_report(
        fs::dir_ls(
          fs::path(data, "primates", species),
          glob = "*_pri_assembly_report.txt"
        )
      )
      # Create a map to rename the NCBI accessions to chromosome names
      refseq_to_chr <- setNames(
        assembly_report$chromosome,
        assembly_report$refseq_accn
      )

      fimo_gr <- load_fimo_results(primate_dir)

      # Rename NCBI accessions
      seqlevels(fimo_gr) <- refseq_to_chr[seqlevels(fimo_gr)]

      fimo_grs[[species]] <- fimo_gr
    }
  }
}

# Function to load the CenSat files
import_censat_bed <- function(path) {
  censat <- rtracklayer::import(path, format = "BED")

  # Drop unnecessary columns
  GenomicRanges::mcols(censat)[, c("score", "itemRgb", "thick")] <- NULL

  # Add region length, will be used for density calculation
  censat$length <- IRanges::width(censat)

  censat
}

load_chm13_censat <- function(chm13_censat_path) {
  chm13_censat <- import_censat_bed(chm13_censat_path)

  # Extract CenSat type
  chm13_censat$type <- tolower(sub("_.*", "", chm13_censat$name))

  # Extract active alpha satellite arrays (end in capital L)
  chm13_censat$type[
    chm13_censat$type == "hor" & grepl("L\\)$", chm13_censat$name)
  ] <- "hor_l"

  chm13_censat
}

load_primate_censat <- function(
  assembly_report_path,
  censat_path,
  censat_mapping
) {
  assembly_report <- load_assembly_report(assembly_report_path)

  primate_censat <- import_censat_bed(censat_path)

  # The primate CenSat files contain also data for unplaced contigs
  # We only keep the canonical chromosomes
  # (those present in the assembly report)
  primate_censat <- primate_censat[
    GenomicRanges::seqnames(primate_censat) %in% assembly_report$chromosome
  ]

  # Extract CenSat type
  primate_censat$type <- tolower(sub("\\(.*", "", primate_censat$name))

  # Keep only CenSat types present in both primates and human assemblies
  primate_censat <- primate_censat[
    primate_censat$type %in% names(censat_mapping)
  ]

  primate_censat
}

# A mapping that goes from abbreviations to full names (used for plotting)
censat_mapping <- c(
  "ct" = "Centromeric transitions",
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
  "active_hor" = "αSat active HOR",
  "rdna" = "rDNA array"
)

# Create a named list to store the CenSat GRanges objects
censat_grs <- list()

for (directory in dir_ls(data, type = "directory")) {
  if (path_file(directory) == "chm13") {
    censat <- dir_ls(directory, glob = "*_censat_v2.1.bed")
    censat_gr <- load_chm13_censat(censat)

    species <- path_file(directory)
    censat_grs[[species]] <- censat_gr

  } else if (path_file(directory) == "primates") {

    for (primate_dir in dir_ls(directory, type = "directory")) {
      report <- dir_ls(primate_dir, glob = "*_pri_assembly_report.txt")
      censat <- dir_ls(primate_dir, glob = "*_CenSat_v2.0.bed")

      censat_gr <- load_primate_censat(report, censat, censat_mapping)

      species <- path_file(primate_dir)
      censat_grs[[species]] <- censat_gr
    }
  }
}

# Create a map to convert GenBank accession names to species names
species_map <- c(
  "PonAbe" = "Pongo abelii",
  "PonPyg" = "Pongo pygmaeus",
  "GorGor" = "Gorilla gorilla gorilla",
  "PanPan" = "Pan paniscus",
  "PanTro" = "Pan troglodytes",
  "SymSyn" = "Symphalangus syndactylus",
  "chm13" = "Homo sapiens"
)

# Create a pre-allocated list to store CTCF density calculation results
density_list <- vector("list", length(names(fimo_grs)))
names(density_list) <- names(fimo_grs)

for (species in names(fimo_grs)) {
  res_gr <- fimo_grs[[species]]
  cen_gr <- censat_grs[[species]]

  mcols(cen_gr)[["overlaps"]] <- countOverlaps(cen_gr, res_gr)
  cen_df <- as.data.frame(mcols(cen_gr))
  cen_df$type <- censat_mapping[cen_df$type]

  density_list[[species]] <- cen_df |>
    dplyr::group_by(type) |>
    dplyr::summarise(
      total_length = sum(length),
      total_overlaps = sum(overlaps)
    ) |>
    dplyr::mutate(density = total_overlaps / (total_length / 1e6)) |>
    dplyr::select(c(type, density)) |>
    tidyr::pivot_wider(names_from = type, values_from = density) |>
    tibble::add_column(species = species_map[[species]], .before = 1)
}

density_df <- as.data.frame(dplyr::bind_rows(density_list))
density_df[density_df == 0] <- NA

# Set species as rownames to be used by ggtree
rownames(density_df) <- density_df$species
density_df <- density_df |> select(-species)

# Define a CenSat order for plotting and reorder the dataframe
censat_order <- c(
  "αSat active HOR", "αSat inactive HOR", "αSat monomeric",
  "αSat divergent", "βSat", "γSat", "HSat1A", "HSat1B",
  "HSat2", "HSat3", "Centromeric transitions",
  "Other satellites", "rDNA array"
)

density_df <- density_df[, censat_order]

phylogenetic_tree_plot <- ggtree(
  tree,
  aes(x, y),
  linewidth = 0.5
) +
  geom_tiplab(
    align = TRUE,
    color = "black",
    size = 2.5
  ) +
  geom_treescale(
    x = 0.0025,
    y = 5,
    fontsize = 2,
    linesize = 0.5
  )
suppressMessages({
  phylogenetic_tree_plot <- gheatmap(
    phylogenetic_tree_plot,
    density_df,
    offset = 0.05,
    width = 6,
    font.size = 2.5,
    colnames_angle = 25,
    hjust = 1,
    color = "black"
  ) +
    scale_fill_viridis_c(
      option = "G",
      trans = "log10",
      na.value = "grey95",
      direction = -1,
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000")
    ) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = 0, r = 0, b = 27, l = 0, unit = "pt")
    ) +
    coord_cartesian(clip = "off")
})

ggsave(
  fs::path(plots, "density_across_primates.svg"),
  plot = phylogenetic_tree_plot,
  width = 6,
  height = 2)