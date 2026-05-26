# Load necessary libraries
suppressPackageStartupMessages({
    suppressWarnings({
        library("rtracklayer", quietly = TRUE, verbose = FALSE)
        library("purrr", quietly = TRUE, verbose = FALSE)
        library("dplyr", quietly = TRUE, verbose = FALSE)
        library("GenomicRanges", quietly = TRUE, verbose = FALSE)
        library("readr", quietly = TRUE, verbose = FALSE)
        library("fs", quietly = TRUE, verbose = FALSE)
        library("yaml", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../../config.yml")

res <- path(config$res) # Path storing general results
fimo_res <- path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data <- path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
plots <- path(config$plots) # Path to store any generated plots

load_fimo_results <- function(genome, fimo_res_dir) {
  # Load all FIMO results for the specified genome
  fimo_files <- dir_ls(path(fimo_res_dir, genome), glob = "*.tsv")
  
  ctcf_dfs <- map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
  ) %>%
    select(-c(motif_alt_id, score, `p-value`, `q-value`, matched_sequence)) # Remove unnecessary columns
  )
  
  res <- bind_rows(ctcf_dfs) # Create one complete dataframe with the results
  
  # Convert the dataframe into a GRanges object for downstream interval manipulation
  ctcf_gr <- makeGRangesFromDataFrame(
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

hg38 <- load_fimo_results("hg38", fimo_res)
chm13 <- load_fimo_results("chm13", fimo_res)

# Convert RefSeq names to UCSC style chromosome naming
hg38_refseq_map <- read_csv(path(data, "hg38", "refseq_to_chr.csv"), show_col_types = FALSE)
chr_map <- setNames(hg38_refseq_map$ucsc, hg38_refseq_map$refseq)

hg38$grl <- endoapply(hg38$grl, function(gr) {
  seqlevels(gr) <- chr_map[seqlevels(gr)]
  gr
})

# Add unique region names for tracking through liftover
hg38$grl <- endoapply(hg38$grl, function(gr) {
  gr$name <- paste0("hg38_region_", seq_along(gr))
  gr
})

chm13$grl <- endoapply(chm13$grl, function(gr) {
  gr$name <- paste0("chm13_region_", seq_along(gr))
  gr
})

# Save the CHM13 and hg38 GRanges objects for each motif to be used for orthology mapping
export_per_motif <- function(grl, assembly, out_dir) {
  for (motif in names(grl)) {
    export(grl[[motif]], path(out_dir, paste0(assembly, "_", motif, ".bed")), format = "BED")
  }
}

export_per_motif(hg38$grl, "hg38", path(res, "misc"))
export_per_motif(chm13$grl, "chm13", path(res, "misc"))