# Load necessary libraries
suppressPackageStartupMessages({
    suppressWarnings({
        library("purrr", quietly = TRUE, verbose = FALSE)
        library("dplyr", quietly = TRUE, verbose = FALSE)
        library("GenomicRanges", quietly = TRUE, verbose = FALSE)
        library("readr", quietly = TRUE, verbose = FALSE)
        library("fs", quietly = TRUE, verbose = FALSE)
        library("yaml", quietly = TRUE, verbose = FALSE)
        library("karyoploteR", quietly = TRUE, verbose = FALSE)
        library("rtracklayer", quietly = TRUE, verbose = FALSE)
    })
})

# Read the configuration YAML file that stores the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

fimo_res = fs::path(config$fimo_res) # Path containing FIMO results for CHM13, HG002, T2T primates, and human pangenome asseblies
data = fs::path(config$data) # Path containing metadata for CHM13, HG002, T2T primates, and human pangenome asseblies
plots = fs::path(config$plots) # Path to store any generated plots

load_fimo_results <- function(genome, fimo_res_dir) {
  # Load all FIMO results for the specified genome
  fimo_files <- dir_ls(fs::path(fimo_res_dir, genome), glob = "*.tsv")
  
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

chm13 <- load_fimo_results("chm13", fimo_res)

# Load CHM13 chromosome sizes 
chm13.chrom.sizes <- read_tsv(
    fs::path(data, "chm13", "chm13v2.0.fa.chrom.sizes"),
    show_col_types = FALSE,
    col_names = c("chr", "size")
)

# Add chromosome start (at 0) for each chromosome
chm13.chrom.sizes <- chm13.chrom.sizes %>%
    mutate(start = 0, .after = chr)

# Load the provided cytobands
chm13.cytobands <- read_tsv(
    fs::path(data, "chm13", "chm13v2.0_cytobands_allchrs.bed"),
    col_names = c("chr", "start", "end", "name", "gieStain"),
    show_col_types = FALSE
)

# Make GRanges objects for both the genome coordinates and the cytobands for karyotype plots 
chm13.genome.gr <- makeGRangesFromDataFrame(
  chm13.chrom.sizes,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "size",
  keep.extra.columns = TRUE
)

chm13.cytobands.gr <- makeGRangesFromDataFrame(
  chm13.cytobands,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE
)

# Create a karyotype plot showing the distribution of CTCF sites across the genome
suppressMessages({
    pdf(fs::path(plots, "karyotype_plot.pdf"), width = 4, height = 11.3)
    
    # Load genome using the T2T-CHM13 genome data
    kp <- plotKaryotype(
        genome = chm13.genome.gr, 
        cytobands = chm13.cytobands.gr, 
        plot.type = 1, # Karyotype with room above for the plot
        cex = 0.85
    )

    # Plot CTCF density per 1Mb
    kp <- kpPlotDensity(
        kp, 
        chm13$gr,
        window.size = 1e6,
        col="#6C6E72"
    )
    dev.off()
})