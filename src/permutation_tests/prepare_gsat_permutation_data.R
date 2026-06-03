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

# Create a dir to store the data that will be used in GAT
gat_dir <- fs::path(res_path, "misc", "gat")
fs::dir_create(gat_dir)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"), format = "BED")

# Extract CenSat type from the names 
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Convert the GRanges object to a dataframe
censat_df <- as.data.frame(chm13_censat_gr)
gsat_df <- censat_df |> filter(type == "gsat") # Keep only gamma-satellites
gsat_df$family <- str_extract(gsat_df$name, "(?<=\\().*?(?=\\))") # Extract gamma-satellite family

# A mapping that normalizes the names of the different gamma-satellite types
gsat_mapping <- c(
  "GSATII" = "GSATII",
  "GSATII/TAR" = "GSATII-TAR composite",
  "GSATX;GSATII" = "GSATII-GSATX composite",
  "GSAT" = "GSAT",
  "GSATII,TAR" = "GSATII-TAR composite",
  "TAR/GSATII" = "GSATII-TAR composite",
  "GSATX" = "GSATX",
  "GSATII,rnd-6_family-1431" = "GSATII",
  "GSATII/TAR1" = "GSATII-TAR composite"
)

# Rename the gamma-satellite families for plotting
gsat_df$family <- gsat_mapping[gsat_df$family]
gsat_df <- gsat_df |> select(seqnames, start, end, family)

# Save the gamma-satellite coordinates as a BED file
gsat_df |>
  write_tsv(
    fs::path(gat_dir, "gsat.bed"),
    col_names = FALSE
  )