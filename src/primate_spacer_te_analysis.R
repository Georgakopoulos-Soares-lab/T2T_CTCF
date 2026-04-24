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
    library("tidyverse", quietly = TRUE, verbose = FALSE)
    library("colorspace", quietly = TRUE, verbose = FALSE)
    library("patchwork", quietly = TRUE, verbose = FALSE)
    library("rlist", quietly = TRUE, verbose = FALSE)
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

spacer_ctcf_te_overlap <- function(
  species,
  full_name,
  censat_bed,
  rmsk_bed,
  regex
) {

  # Load the FIMO results into a GRanges object
  gr <- load_fimo_results(
    fs::path(fimo_res, "primates", species)
  )

  # Load the CenSat annotations
  censat_gr <- rtracklayer::import(
    fs::path(data, "primates", species, censat_bed),
    format = "BED"
  )

  # Add a column specifing the type of the CenSat region
  censat_gr$type <- tolower(sub(regex, "", censat_gr$name))

  # Drop unused columns and extract the
  # subterminal heterochromatin (subterm) regions
  GenomicRanges::mcols(censat_gr)[, c("score", "itemRgb", "thick")] <- NULL
  subterm_gr <- censat_gr[censat_gr$type == "subterm"]

  # Create a broader subterm GRanges object by combining regions 1Mb apart
  # This is done so that we can identify the spacers between regions
  combined_subterm_gr <- GenomicRanges::reduce(
    subterm_gr,
    min.gapwidth = 1000001
  )

  # Spacers are found by removing the original
  # subterm regions from the combined ones
  spacer_gr <- GenomicRanges::setdiff(combined_subterm_gr, subterm_gr)

  # Read the RepeatMasker data
  rmsk_df <- readr::read_tsv(
    fs::path(
      data,
      "primates",
      species,
      rmsk_bed
    ),
    col_names = c(
      "seqnames", "start", "end",
      "name", "score", "strand",
      "thickStart", "thickEnd",
      "itemRgb", "blockCount",
      "blockSizes", "blockStarts",
      "RM index", "orig_out"
    ),
    show_col_types = FALSE
  )

  # Select only the relevant columns
  rmsk_df <- rmsk_df |>
    dplyr::select(
      c(
        "seqnames", "start", "end",
        "name", "score", "strand"
      )
    )

  # Split the combined names to name/class
  rmsk_df <- rmsk_df |>
    tidyr::separate_wider_delim(
      name,
      delim = "#",
      names = c("name", "class"),
      too_few = "align_start"
    )

  # Split the combined class to class/subclass
  rmsk_df <- rmsk_df |>
    tidyr::separate_wider_delim(
      class,
      delim = "/",
      names = c("class", "subclass"),
      too_few = "align_start"
    )

  # Convert the RepeatMasker dataframe to a GRanges object
  rmsk_gr <- GenomicRanges::makeGRangesFromDataFrame(
    rmsk_df,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )

  # Find which CTCF sites fall in spacer ranges
  spacer_hits <- GenomicRanges::findOverlaps(
    gr,
    spacer_gr,
    ignore.strand = TRUE
  )
  spacer_ctcf_sites <- GenomicRanges::sort(
    gr[S4Vectors::queryHits(spacer_hits)]
  )

  # Find which of the spacer CTCF sites overlap with RepeatMasker regions
  spacer_rmsk_hits <- GenomicRanges::findOverlaps(
    spacer_ctcf_sites,
    rmsk_gr,
    ignore.strand = TRUE
  )

  spacer_ctcf_rmsk <- spacer_ctcf_sites[
    S4Vectors::queryHits(spacer_rmsk_hits)
  ]

  # Add the RepeatMasker class, subclass, and name to each CTCF site
  spacer_ctcf_rmsk$class  <- rmsk_gr$class[
    S4Vectors::subjectHits(spacer_rmsk_hits)
  ]

  spacer_ctcf_rmsk$subclass <- rmsk_gr$subclass[
    S4Vectors::subjectHits(spacer_rmsk_hits)
  ]

  spacer_ctcf_rmsk$name <- rmsk_gr$name[
    S4Vectors::subjectHits(spacer_rmsk_hits)
  ]

  # Convert the GRanges object to a dataframe and add the corresponding species
  spacer_ctcf_rmsk_df <- as.data.frame(spacer_ctcf_rmsk)

  spacer_ctcf_rmsk_df <- spacer_ctcf_rmsk_df |>
    dplyr::mutate(
      species = full_name
    )

  list(
    total_spacer_sites = length(spacer_ctcf_sites),
    rmsk_spacer_sites = spacer_ctcf_rmsk_df
  )
}

# Detect overlaps between subterminal spacer CTCF sites
# and transposable elements
suppressWarnings({
  gorgor_spacer_overlaps <- spacer_ctcf_te_overlap(
    "GorGor",
    "Gorilla gorilla gorilla",
    "mGorGor1_v2.0_CenSat_v2.0.bed",
    "mGorGor1_v2.0_pri.RepeatMasker_v1.1.bed",
    "\\(.*"
  )

  panpan_spacer_overlaps <- spacer_ctcf_te_overlap(
    "PanPan",
    "Pan paniscus",
    "mPanPan1_v2.0_CenSat_v2.0.bed",
    "mPanPan1_v2.0_pri.RepeatMasker_v1.1.bed",
    "\\(.*"
  )

  pantro_spacer_overlaps <- spacer_ctcf_te_overlap(
    "PanTro",
    "Pan troglodytes",
    "mPanTro3_v2.0_CenSat_v2.0.bed",
    "mPanTro3_v2.0_pri.RepeatMasker_v1.1.bed",
    "\\(.*"
  )

  symsyn_spacer_overlaps <- spacer_ctcf_te_overlap(
    "SymSyn",
    "Symphalangus syndactylus",
    "mSymSyn1_v2.1_CenSat_v2.0.bed",
    "mSymSyn1_v2.0_pri.RepeatMasker_v1.1.bed",
    "_.*"
  )
})

# Create a list with the spacer-transposable element ovelaps
overlaps <- list(
  gorgor_spacer_overlaps,
  panpan_spacer_overlaps,
  pantro_spacer_overlaps,
  symsyn_spacer_overlaps
)

# Here, we create a dataframe containing the percentage of spacer CTCF
# sites which overlap with transposable elements for each species
spacer_ctcf_te_perc <- map_dfr(overlaps, \(obj) {
  # Count unique sites in transposable elements
  # (we want unique because GRanges::findOverlaps() allows duplicate overlaps)
  in_te <- nrow(
    distinct(
      obj$rmsk_spacer_sites, seqnames, start, end, width, strand, motif_id
    )
  ) /
    obj$total_spacer_sites * 100
  tibble(
    species = unique(obj$rmsk_spacer_sites$species),
    in_te = in_te,
    not_in_te = 100 - in_te
  )
}) |>
  pivot_longer(
    c(in_te, not_in_te),
    names_to = "category",
    values_to = "percent"
  ) |>
  mutate(
    category = factor(category, levels = c("not_in_te", "in_te"))
  )

# Here we create a stacked bar chart showing the above percentages
# and save to a PDF file
stacked_perc <- ggplot(
  spacer_ctcf_te_perc,
  aes(x = species, y = percent, fill = category)
) +
  geom_col(
    color = "black",
    linewidth = 0.5,
    width = 0.7
  ) +
  labs(
    x = NULL,
    y = "Percentage (%)",
    fill = "Category"
  ) +
  scale_x_discrete(
    labels = scales::label_wrap(20),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_fill_manual(
    values = c(not_in_te = "grey70", in_te = "#9D5237")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      face = "italic",
      size = 8
    ),
    axis.text.y = element_text(
      size = 6
    ),
    axis.title.y = element_text(
      size = 8
    ),
    legend.position = "none",
    plot.background = NULL,
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
  )

ggsave(
  fs::path(plots, "primate_spacer_stacked.pdf"),
  stacked_perc,
  width = 3,
  height = 1.5,
  bg = "transparent"
)

# Now we compute transposable element distribution for the CTCF spacer sites
spacer_ctcf_te_dist <- map_dfr(overlaps, \(obj) {
  as.data.frame(
    (table(obj$rmsk_spacer_sites$subclass) / nrow(obj$rmsk_spacer_sites)) * 100
  ) |>
    dplyr::rename(
      transposable_element = Var1,
      percentage = Freq
    ) |>
    dplyr::mutate(
      species = unique(obj$rmsk_spacer_sites$species)
    ) |>
    # Attach the parent class for each subclass
    dplyr::left_join(
      obj$rmsk_spacer_sites |>
        dplyr::select(c(class, subclass)) |>
        distinct(),
      by = c("transposable_element" = "subclass")
    )
}) |>
  # Ensure every species x TE combination exists; missing ones get 0%
  complete(
    transposable_element,
    species,
    fill = list(percentage = 0)
  ) |>
  # Propagate class to rows added by complete() (which lack class info)
  dplyr::group_by(transposable_element) |>
  tidyr::fill(class, .direction = "downup") |>
  dplyr::ungroup()

# Order TEs by class, then rename satellite subclasses for display
spacer_ctcf_te_dist <- spacer_ctcf_te_dist |>
  dplyr::arrange(class, transposable_element) |>
  dplyr::mutate(
    transposable_element = case_match(
      transposable_element,
      "acro" ~ "Acrocentric",
      "centr" ~ "(Peri)centromeric",
      "subtelo" ~ "Subtelomeric",
      .default = transposable_element
    )
  ) |>
  # Lock factor order (reversed so first class appears at top of plot)
  dplyr::mutate(
    transposable_element = forcats::fct_inorder(transposable_element) |>
      forcats::fct_rev()
  )

# Compute y-axis extents per class for drawing
# the right-side bracket annotations
class_ranges <- spacer_ctcf_te_dist |>
  dplyr::distinct(transposable_element, class) |>
  dplyr::mutate(y = as.integer(transposable_element)) |>
  dplyr::summarise(
    ymin = min(y) - 0.35,
    ymax = max(y) + 0.35,
    ymid = (min(y) + max(y)) / 2,
    .by = class
  )

# X-coordinates for bracket geometry and class labels (outside the plot area)
x_max <- 98
x_bracket <- x_max * 1.065
x_tick <- x_max * 0.015
x_label <- x_max * 1.08

# Here, we define a palette for the species
species_pal <- c(
  "Gorilla gorilla gorilla" = "#4E7B8A",
  "Pan paniscus" = "#C98A5F",
  "Pan troglodytes" = "#8A9B6E",
  "Symphalangus syndactylus" = "#8B5E83"
)

# We create a grouped bar chart showing the percentage of CTCF spacer
# site assigned to each transposable element subclass
grouped_te_barplot <- ggplot(
  spacer_ctcf_te_dist,
  aes(y = transposable_element, x = percentage, fill = species)
) +
  geom_col(position = "dodge", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = species_pal) +
  geom_segment(
    data = class_ranges,
    aes(x = x_bracket, xend = x_bracket, y = ymin, yend = ymax),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_segment(
    data = class_ranges,
    aes(x = x_bracket - x_tick, xend = x_bracket + 0.38, y = ymin, yend = ymin),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_segment(
    data = class_ranges,
    aes(x = x_bracket - x_tick, xend = x_bracket + 0.38, y = ymax, yend = ymax),
    inherit.aes = FALSE, linewidth = 0.4
  ) +
  geom_text(
    data = class_ranges,
    aes(x = x_label, y = ymid, label = class),
    inherit.aes = FALSE, hjust = 0, size = 2
  ) +
  labs(y = "Transposable elements", x = "Percentage (%)", fill = "Species") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position    = "none",
    panel.grid.major.x = element_line(color = "gray70", linewidth = 0.3),
    plot.background    = NULL,
    axis.text.y        = element_text(size = 6),
    axis.text.x        = element_text(size = 6),
    axis.title.x       = element_text(size = 8),
    axis.title.y       = element_text(size = 8),
    plot.margin        = margin(5, 28, 5, 5)
  )

ggsave(
  fs::path(plots, "primate_spacer_te_grouped.pdf"),
  grouped_te_barplot,
  width = 3,
  height = 3.3,
  bg = "transparent"
)