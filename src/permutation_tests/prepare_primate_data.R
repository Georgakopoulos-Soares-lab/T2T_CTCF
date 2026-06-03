# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
    library("plyranges", quietly = TRUE, verbose = FALSE)
    library("Biostrings", quietly = TRUE, verbose = FALSE)
    library("Rsamtools", quietly = TRUE, verbose = FALSE)
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

# Create a dir to store the data that will be used in GAT
primates_gat_dir <- fs::path(res_path, "misc", "gat", "primates")
fs::dir_create(primates_gat_dir)

# Helper function to load FIMO results
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

  # Remove strand information
  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
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

# Helper function to export a BED file for the GAT permutation analysis
export_gat_bed <- function(
  gr,
  species,
  filename
) {
  # Define the directory to save the FIMO BED file
  out_dir <- fs::path(primates_gat_dir, species)

  # Create the directory if it does not exist
  fs::dir_create(out_dir)

  # Export the BED file
  rtracklayer::export(
    gr,
    fs::path(out_dir, filename),
    format = "BED"
  )
}

# Export the CTCF BED file containing FIMO results for each species
purrr::iwalk(
  fimo_grs, \(gr, s) export_gat_bed(gr, s, "ctcf_sites.bed")
)

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
  chm13_censat_gr |>
    plyranges::select(name)
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
      sub("subTerm_mon*", "subterm", primate_censat$name)
    )
  )

  # Keep only regions of interest
  primate_censat <- primate_censat[
    primate_censat$name %in% names(censat_mapping)
  ]

  # Remove the unnecessary metadata columns for the GAT analysis
  primate_censat |>
    plyranges::select(name)
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

# Export the formatted CenSat BED file for each species
purrr::iwalk(
  censat_grs, \(gr, s) export_gat_bed(gr, s, "censat.bed")
)

# Helper function which defines the workspace and isochores
# across primate assemblies to be used in the permutation based testing
create_permutation_data <- function(
  chrom_sizes_path,
  fasta_path
) {
  # Load the chromosome sizes which will be used a the
  # workspace in the permutation analysis
  workspace <- readr::read_tsv(
    chrom_sizes_path,
    show_col_types = FALSE,
    col_names = c("chrom", "size")
  ) |>
    dplyr::mutate(start = 0, workspace = "ws") |>
    dplyr::filter(!chrom %in% c("chrM", "MT")) |>
    dplyr::select(chrom, start, size, workspace)

  # Create a named vector that will be used to tile the genome
  lengths <- setNames(
    workspace$size,
    workspace$chrom
  )

  # Tile the genome in 100 Kb windows
  windows <- GenomicRanges::tileGenome(
    lengths,
    tilewidth = 100e3,
    cut.last.tile.in.chrom = TRUE
  )

  # Provide the FASTA file as a pointer for R to use
  # Make sure the genome file is indexed using Samtools
  fa <- Rsamtools::FaFile(fasta_path)

  # Get the underlying sequence for each window
  window_sequences <- Biostrings::getSeq(fa, windows)

  # Calculate GC content (counts / region's length)
  # and add it as metadata to the windows GRanges
  gc_count <- as.integer(
    Biostrings::letterFrequency(
      window_sequences, "GC"
    )
  )

  # Create a vector for tile widths
  w <- GenomicRanges::width(windows)

  # Add GC count and GC content to each tile
  windows <- windows |>
    plyranges::mutate(
      gc_count = gc_count,
      gc_content = gc_count / w
    )

  # Define isochores to be used with GAT
  isochores <- windows |>
    plyranges::mutate(
      name = cut(
        gc_content,
        breaks = quantile(gc_content, probs = seq(0, 1, 1 / 6), na.rm = TRUE),
        labels = c("L1", "L2", "M1", "M2", "H1", "H2"),
        include.lowest = TRUE
      )
    ) |>
    plyranges::select(-c(gc_content, gc_count))

  # Return both the workspace and the isochores
  list(workspace = workspace, isochores = isochores)
}

# Create tiles for the CHM13 genome
chm13_tiles <- create_permutation_data(
  fs::path(data_path, "chm13", "chm13v2.0.chrom.sizes"),
  fs::path(data_path, "chm13", "chm13v2.0.fa.gz")
)

# Update the primate paths with the paths for the chromosome
# sizes and the genome FASTA files
primate_paths <- primate_paths |>
  dplyr::mutate(
    chrom_sizes = purrr::map_chr(
      primate_dirs,
      \(d) fs::dir_ls(d, glob = "*.chrom.sizes")
    ),
    fasta = purrr::map_chr(
      primate_dirs,
      \(d) fs::dir_ls(d, glob = "*.fna.gz")
    )
  )

# Create tiles for the primates genomes and also create a list
# containing all results
tiles <- c(
  list(chm13 = chm13_tiles),
  purrr::map2(
    primate_paths$chrom_sizes,
    primate_paths$fasta,
    create_permutation_data
  ) |>
    setNames(tolower(primate_paths$species))
)

# Save the GAT workspace and the isochore regions for GAT analysis
purrr::iwalk(tiles, \(tile, species) {
  out_dir <- fs::path(primates_gat_dir, species)
  fs::dir_create(out_dir)

  readr::write_tsv(
    tile$workspace,
    fs::path(out_dir, "workspace.bed"),
    col_names = FALSE
  )

  export_gat_bed(tile$isochores, species, "isochores.bed")
})