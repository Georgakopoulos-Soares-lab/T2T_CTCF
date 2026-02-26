# Load necessary libraries
suppressPackageStartupMessages({
    suppressWarnings({
        library("purrr", quietly = TRUE, verbose = FALSE)
        library("regioneR", quietly = TRUE, verbose = FALSE)
        library("dplyr", quietly = TRUE, verbose = FALSE)
        library("GenomicRanges", quietly = TRUE, verbose = FALSE)
        library("readr", quietly = TRUE, verbose = FALSE)
        library("fs", quietly = TRUE, verbose = FALSE)
        library("yaml", quietly = TRUE, verbose = FALSE)
        library("rtracklayer", quietly = TRUE, verbose = FALSE)
        library("tidyr", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

fimo_res <- fs::path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data <- fs::path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
res_path <- fs::path(config$res) # Path where all results are stored

load_fimo_results <- function(genome, fimo_res_dir) {
  # Load all FIMO results for the specified genome
  fimo_files <- fs::dir_ls(fs::path(fimo_res_dir, genome), glob = "*.tsv")

  ctcf_dfs <- purrr::map(fimo_files, ~ read_tsv(
    .x,
    comment = "#", # Remove comment lines at the end of each FIMO file
    show_col_types = FALSE
  ) |>
    select(-c(motif_alt_id, score, `p-value`, `q-value`, matched_sequence)) # Remove unnecessary columns
  )

  res <- dplyr::bind_rows(ctcf_dfs) # Create one complete dataframe with the results

  # Convert the dataframe into a GRanges object for downstream interval manipulation
  ctcf_gr <- GenomicRanges::makeGRangesFromDataFrame(
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

chm13 <- load_fimo_results("chm13", fimo_res)

# A mapping that goes from abbreviations to full names (used for plotting)
censat_mapping <- c(
    "ct" = "Centromeric transitions",
    "censat" = "Other satellites",
    "gsat" = "γSat",
    "bsat" = "βSat",
    "hsat3" = "HSat3",
    "hsat2" = "HSat2",
    "hsat1B" = "HSat1B",
    "hsat1A" = "HSat1A",
    "mon" = "αSat monomeric",
    "dhor" = "αSat divergent",
    "hor" = "αSat inactive HOR",
    "hor_l" = "αSat active HOR",
    "rDNA" = "rDNA array"
)

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"), format = "BED")

# Extract CenSat type from the names 
chm13_censat_gr$type <- sub("_.*", "", chm13_censat_gr$name)

# If an HOR has L at the end of the name then this indicates that this is the active alpha-satellite HOR
chm13_censat_gr$type[chm13_censat_gr$type == "hor" & grepl("L\\)$", chm13_censat_gr$name)] <- "hor_l"

# Drop unnecessary collumns
mcols(chm13_censat_gr)[, c("score", "itemRgb", "thick")] <- NULL

# Read the chromosome sizes for the CHM13 genome
chm13_chrom_sizes <- read_tsv(
    fs::path(data, "chm13", "chm13v2.0.fa.chrom.sizes"), 
    show_col_types = FALSE,
    col_names = c("chrom", "size"))

# Create the genome object used in regioneR and store it once
# This is memoised and then loaded intantaneously across the many runs
genome_and_mask <- getGenomeAndMask(as.data.frame(chm13_chrom_sizes), mask = NA)

# We define a function to run the permutation testing we perform to test if the 
# enrichment of CTCF sites within certain regions is statistically significant
# The core idea is that we take the original FIMO CTCF hits per motif and randomly
# distribute them across the genome. We restrict the random shuffling to be per chromosome
run_permtest <- function(motif, region) {
  
  # Remove the strand information so that we perform a strand-agnostic overlap
  strand(motif) <- "*"
  strand(region) <- "*"
  
  pt <- permTest(
    A = motif, # Motif FIMO res
    B = region, # Regions assigned to one CenSat type
    ntimes = 10000, # Number of permutations
    randomize.function = randomizeRegions,
    evaluate.function = numOverlaps,
    alternative = "greater", # We only test for enrichment (not depletion)
    genome = genome_and_mask$genome,
    mask = NULL,
    allow.overlaps = TRUE,
    per.chromosome = TRUE,
    verbose = FALSE,
    mc.set.seed = FALSE, # We remove the seed here, so that each core uses the same seed defined in the start
    force.parallel = TRUE,
    mc.cores = 130
  )
  
  # Store the results in a dataframe 
  data.frame(
    pval     = pt$numOverlaps$pval,
    zscore   = pt$numOverlaps$zscore
  )
}

set.seed(42) # Set seed so that each permutation test result is reproducible
results_list <- list() # Define a list where we will store the results

# Run the permutation testing for each motif and for each CenSat type.
# We store the results for each combination in the list
for (type in unique(chm13_censat_gr$type)) {
  censat_type <- chm13_censat_gr[chm13_censat_gr$type == type]
  for (motif_name in names(chm13$grl)) {
    message(sprintf("Running: %s x %s", type, motif_name))
    motif_gr <- chm13$grl[[motif_name]]
    res <- run_permtest(motif_gr, censat_type)
    res$type  <- type
    res$motif <- motif_name
    results_list[[length(results_list) + 1]] <- res
  }
}

results <- do.call(rbind, results_list)
results$padj <- p.adjust(results$pval, method = "BH")
# Save results to a CSV file 
# We will load these results in the code we use for the plotting of the densities across CenSat types
write.csv(results, fs::path(res_path, "misc", "permtest_results.csv"), row.names = FALSE)