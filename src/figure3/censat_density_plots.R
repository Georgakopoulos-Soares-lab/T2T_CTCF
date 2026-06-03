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

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res_path <- fs::path(config$res) # Path containing all results
# Path containing FIMO results for T2T primates
fimo_res <- fs::path(res_path, "fimo_res")

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

# Directory containing the GAT permutation results
primates_gat_dir <- fs::path(res_path, "misc", "gat", "primates")

# Load primate & human phylogenetic tree
tree <- ape::read.tree(
  fs::path(
    data_path,
    "multiple_sequence_alignments",
    "8-t2t-apes-2023v2.nh"
  )
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

# Helper function to load FIMO results
load_fimo_results <- function(genome_res) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(genome_res, glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "# ", # Remove comment lines at the end of each FIMO file
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

# Helper function that loads the assembly report for each primate species
# to find which chromosomes are the canonical ones
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

# Loop over the different directories in the FIMO results path
# We need only the CHM13 and primates folders
target_dirs <- c(
  fs::dir_ls(
    fimo_res,
    regexp = "chm13$",
    type = "directory"
  ),
  fs::dir_ls(
    fimo_res / "primates",
    type = "directory"
  )
)

# Load results, naming each element by its species folder
fimo_grs <- target_dirs |>
  purrr::set_names(~ tolower(fs::path_file(.x))) |>
  purrr::map(load_fimo_results)

# Helper function to load the CenSat data for CHM13
load_chm13_censat <- function(chm13_censat_path) {

  # Read the CenSat BED file into a GRanges object
  chm13_censat_gr <- rtracklayer::import(chm13_censat_path, format = "BED")

  # Extract CenSat type
  chm13_censat_gr$name <- tolower(sub("_.*", "", chm13_censat_gr$name))

  # Define active alpha satellites based on CenSat color
  chm13_censat_gr$name[
    chm13_censat_gr$name == "hor" & chm13_censat_gr$itemRgb == "#990000"
  ] <- "hor_l"

  # Remove the unnecessary metadata columns for the GAT analysis
  chm13_censat_gr <- chm13_censat_gr |>
    plyranges::select(name)

  chm13_censat_gr$name <- censat_mapping[chm13_censat_gr$name]
  chm13_censat_gr$length <- GenomicRanges::width(chm13_censat_gr)

  chm13_censat_gr
}

# Helper function to load the CenSat data for primate species
load_primate_censat <- function(
  assembly_report_path,
  censat_path,
  censat_mapping
) {

  # Load the NCBI assembly report
  assembly_report <- load_assembly_report(assembly_report_path)

  # Read the CenSat BED file into a GRanges object
  primate_censat <- rtracklayer::import(censat_path, format = "BED")

  # The primate CenSat files contain also data for unplaced contigs
  # We only keep the canonical chromosomes
  # (those present in the assembly report)
  primate_censat <- primate_censat[
    GenomicRanges::seqnames(primate_censat) %in% assembly_report$chromosome
  ]

  # Extract CenSat type
  primate_censat$name <- tolower(
    sub("\\(.*", "",
      sub("subTerm_mon.*", "subterm", primate_censat$name)
    )
  )

  # Keep only regions of interest
  primate_censat <- primate_censat[
    primate_censat$name %in% names(censat_mapping)
  ]

  # Remove the unnecessary metadata columns for the GAT analysis
  primate_censat <- primate_censat |>
    plyranges::select(name)

  primate_censat$name <- censat_mapping[primate_censat$name]
  primate_censat$length <- GenomicRanges::width(primate_censat)

  primate_censat
}

# Define the directories that contain data for each primate species
primate_dirs <- fs::dir_ls(data_path / "primates", type = "directory")

# Remove the AQERs one, since this contains certain accelerated elements
primate_dirs <- primate_dirs[path_file(primate_dirs) != "AQERs"]

# Populate paths with the specific files needed for the analysis
# We need the assembly reports and the CenSat annotations
primate_paths <- tidyr::tibble(
  species = tolower(fs::path_file(primate_dirs)),
  report = purrr::map_chr(
    primate_dirs,
    \(d) dir_ls(d, glob = "*_pri_assembly_report.txt")
  ),
  censat = purrr::map_chr(
    primate_dirs,
    \(d) dir_ls(d, glob = "*_CenSat_v2.0.bed")
  )
)

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
  "rdna" = "rDNA array",
  "subterm" = "Subterminal heterochromatin"
)

# Create a list with the CenSat GRanges objects for CHM13 and the 6 primates
censat_grs <- c(
  list(
    chm13 = load_chm13_censat(
      fs::path(data_path, "chm13", "chm13v2.0_censat_v2.1.bed")
    )
  ),
  with(
    primate_paths,
    set_names(
      purrr::pmap(
        list(report, censat, list(censat_mapping)),
        load_primate_censat
      ),
      species
    )
  )
)

# Create a map to convert GenBank accession names to species names
species_map <- c(
  "ponabe" = "Pongo abelii",
  "ponpyg" = "Pongo pygmaeus",
  "gorgor" = "Gorilla gorilla gorilla",
  "panpan" = "Pan paniscus",
  "pantro" = "Pan troglodytes",
  "symsyn" = "Symphalangus syndactylus",
  "chm13" = "Homo sapiens"
)

# Create a pre-allocated list to store CTCF density calculation results
density_list <- vector("list", length(names(fimo_grs)))
names(density_list) <- names(fimo_grs)

# Here, we iterate over species extract the relevant FIMO data and
# CenSat annotations, find how many CTCF sites overlap each annotation
# calculate CTCF density and pivot the data to be accepted by ggtree
suppressWarnings({
  for (species in names(fimo_grs)) {
    res_gr <- fimo_grs[[species]]
    cen_gr <- censat_grs[[species]]

    mcols(cen_gr)[["overlaps"]] <- countOverlaps(cen_gr, res_gr)
    cen_df <- as.data.frame(mcols(cen_gr))

    density_list[[species]] <- cen_df |>
      dplyr::group_by(name) |>
      dplyr::summarise(
        total_length = sum(length),
        total_overlaps = sum(overlaps)
      ) |>
      dplyr::mutate(density = total_overlaps / (total_length / 1e6)) |>
      dplyr::select(c(name, density)) |>
      tidyr::pivot_wider(names_from = name, values_from = density) |>
      tibble::add_column(species = species_map[[species]], .before = 1)
  }
})

# Combine all species results into 1 dataframe
density_df <- as.data.frame(dplyr::bind_rows(density_list))

# Set densities of 0 to be NA so that the log-scale works when plotting
density_df[density_df == 0] <- NA

# # Set species as rownames to be used by ggtree
rownames(density_df) <- density_df$species
density_df <- density_df |> select(-species)

# Define a CenSat order for plotting and reorder the dataframe
censat_order <- c(
  "αSat active HOR", "αSat inactive HOR", "αSat monomeric",
  "αSat divergent", "βSat", "γSat", "HSat1A", "HSat1B",
  "HSat2", "HSat3", "Centromeric transitions",
  "Other satellites", "rDNA array", "Subterminal heterochromatin"
)

# Rearrange the density dataframe to follow the CenSat order
density_df <- density_df[, censat_order]

# We load the results of GAT for each species, select only needed columns
# and keep only enrichment associations. We also add star annotations based
# on q-values
gat_res <- purrr::map(names(species_map), \(sp) {
  read_tsv(
    fs::path(primates_gat_dir, sp, "GAT_results.tsv"),
    show_col_types = FALSE
  ) |>
    dplyr::select(track, annotation, fold, qvalue) |>
    dplyr::filter(fold > 1, qvalue < 0.05) |>
    dplyr::mutate(
      species = sp,
      sig = case_when(
        qvalue < 0.001 & fold > 1 ~ "***",
        qvalue < 0.01 & fold > 1 ~ "**",
        qvalue < 0.05 & fold > 1 ~ "*",
        .default = " "
      )
    )
}) |>
  purrr::list_rbind()

# We make sure that every species x CenSat annotation exists and
# fill the missing ones with NA. Also, we format the q-value star
# order so that it goes from MA1930.2 to MA1929.2 and MA0139.2
gat_res <- gat_res |>
  dplyr::mutate(
    motif = factor(
      track, levels = c("MA1930.2", "MA1929.2", "MA0139.2")
    )
  ) |>
  tidyr::complete(
    tidyr::nesting(species, annotation),
    motif,
    fill = list(sig = " ")
  ) |>
  dplyr::arrange(species, annotation, motif) |>
  dplyr::group_by(species, annotation) |>
  dplyr::summarise(sig = paste(sig, collapse = "\n"), .groups = "drop")

# We add the full names for each CenSat annotation and for each species acronym
gat_res$annotation <- censat_mapping[gat_res$annotation]
gat_res$species <- species_map[gat_res$species]

# We create a plot showing a phylogenetic tree with a heatmap
# depicting CTCF density across CenSat types
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
    offset = 0.06,
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
    coord_cartesian(clip = "off")

  # We extract the coordinates for each tile of the heatmap so that we can
  # overlay the q-value stars on top
  is_hm <- vapply(
    phylogenetic_tree_plot$layers,
    \(l) {
      inherits(l$geom, "GeomTile") &&
        all(c("variable", "value") %in% names(l$data))
    },
    logical(1)
  )

  cell_coords <- phylogenetic_tree_plot$layers[[which(is_hm)[1]]]$data |>
    dplyr::distinct(
      species = lab,
      annotation = as.character(variable),
      x, y
    )

  star_df <- dplyr::inner_join(
    cell_coords, gat_res,
    by = c("species", "annotation")
  )

  phylogenetic_tree_plot <- phylogenetic_tree_plot +
    geom_text(
      data = star_df,
      aes(x = x, y = y - 0.09, label = sig),
      inherit.aes = FALSE,
      size       = 1.5,
      lineheight = 0.75,
      color      = "white"
    ) +
    theme(
      legend.position = "none",
      plot.background = NULL,
      plot.margin = margin(t = 0, r = 0, b = 35, l = 0, unit = "pt")
    )
})

# Save the plot
ggsave(
  fs::path(plots, "density_across_primates.svg"),
  plot = phylogenetic_tree_plot,
  width = 5.5,
  height = 2,
  bg = "transparent"
)