# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("ggplot2", quietly = TRUE, verbose = FALSE)
    library("dplyr", quietly = TRUE, verbose = FALSE)
    library("GenomicRanges", quietly = TRUE, verbose = FALSE)
    library("fs", quietly = TRUE, verbose = FALSE)
    library("yaml", quietly = TRUE, verbose = FALSE)
    library("rtracklayer", quietly = TRUE, verbose = FALSE)
    library("readr", quietly = TRUE, verbose = FALSE)
    library("purrr", quietly = TRUE, verbose = FALSE)
    library("circlize", quietly = TRUE, verbose = FALSE)
    library("tidyr", quietly = TRUE, verbose = FALSE)
    library("stringr", quietly = TRUE, verbose = FALSE)
    library("ComplexHeatmap", quietly = TRUE, verbose = FALSE)
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

  GenomicRanges::strand(fimo_gr) <- "*"

  # Return the FIMO results
  fimo_gr
}

# Load the CHM13 FIMO results
chm13_gr <- load_fimo_results(fs::path(fimo_res, "chm13"))

# Load CHM13 CenSat annotations
chm13_censat_gr <- import(
  fs::path(data, "chm13", "chm13v2.0_censat_v2.1.bed"),
  format = "BED"
)

# Extract CenSat type from the names
chm13_censat_gr$type <- tolower(sub("_.*", "", chm13_censat_gr$name))
chm13_censat_gr$annotation <- tolower(stringr::str_extract(
  chm13_censat_gr$name, "(?<=\\().*?(?=\\))"
))

chm13_censat_gr$annotation[is.na(chm13_censat_gr$annotation)] <- "No annotation"

# Define active alpha satellites based on CenSat color
chm13_censat_gr$type[
  chm13_censat_gr$type == "hor" & chm13_censat_gr$itemRgb == "#990000"
] <- "hor_l"

# Extract the acrocentric short arm regions for plotting
acrocentric_regions <- subset(
  chm13_censat_gr,
  seqnames %in% c("chr13", "chr14", "chr15", "chr21", "chr22") & type == "hor_l"
)

# Set all starts to 1
start(acrocentric_regions) <- 1

# Create one unified region for each chromosome
acrocentric_regions <- unlist(
  range(
    split(
      acrocentric_regions,
      seqnames(acrocentric_regions)
    )
  )
)

acro_ctcf <- as.data.frame(subsetByOverlaps(
  chm13_gr,
  acrocentric_regions,
  ignore.strand = TRUE
)) |>
  mutate(
    query_name = paste0("CTCF_", row_number())
  ) |>
  select(c(seqnames, start, end, query_name, strand))

files <- fs::dir_ls(
  fs::path(res_path, "misc"),
  glob = "*acro_short_arms.position.chr*.tsv"
)

odgi_pos_res <- purrr::map(
  files,
  ~ read_tsv(
    .x,
    skip = 1,
    show_col_types = FALSE,
    col_names = c(
      "query_chrom", "query_start", "query_end",
      "query_name", "query_strand",
      "ref_start", "ref_end", "ref_strand"
    )
  )
) |>
  purrr::list_rbind() |>
  dplyr::mutate(
    query_chrom = str_split_i(query_chrom, "#", 1),
    ref_chrom = str_split_i(ref_start, "#", 1),
    ref_start = as.numeric(str_split_i(ref_start, ",", 2)),
    ref_end = as.numeric(str_split_i(ref_end, ",", 2)),
    ref_start   = pmin(ref_start, ref_end),
    ref_end     = pmax(ref_start, ref_end),
    query_width = query_end - query_start,
    ref_width = ref_end - ref_start
  ) |>
  dplyr::arrange(
    query_chrom, query_start, query_end,
    ref_chrom, ref_start, ref_end,
  ) |>
  dplyr::filter(
    query_width == ref_width
  ) |>
  dplyr::select(
    query_chrom, query_start, query_end,
    query_strand, query_width, ref_chrom,
    ref_start, ref_end, ref_strand, ref_width
  )

target_gr <- GRanges(
  seqnames = odgi_pos_res$ref_chrom,
  ranges   = IRanges(
    start = odgi_pos_res$ref_start,
    end = odgi_pos_res$ref_end
  )
)

hits <- findOverlaps(
  chm13_gr,
  target_gr,
  type = "within",
  ignore.strand = TRUE
)

odgi_pos_res$is_ctcf <- as.integer(
  seq_along(target_gr) %in% subjectHits(hits)
)

odgi_pos_res <- odgi_pos_res |>
  filter(is_ctcf == 1)

stats_per_site <- odgi_pos_res |>
  dplyr::group_by(
    query_chrom,
    query_start,
    query_end
  ) |>
  dplyr::summarise(
    n_targets = n(),
    is_ctcf = as.integer(any(is_ctcf == 1)),
    .groups = "drop"
  )

stats_per_site <- stats_per_site |>
  dplyr::mutate(
    category = case_when(
      n_targets == 1 ~ "1 arm",
      TRUE ~ paste0(n_targets, " arms")
    )
  )

summary_df <- stats_per_site |>
  dplyr::count(category) |>
  dplyr::mutate(pct = (n / sum(n)) * 100) |>
  dplyr::mutate(
    color = adjustcolor(
      dplyr::recode(
        category,
        "1 arm"  = "#FCBBA1",
        "2 arms" = "#FB6A4A",
        "3 arms" = "#DE2D26",
        "4 arms" = "#A50F15",
        "5 arms" = "#67000D",
        .default = "grey90"
      ),
      alpha.f = 1
    )
  ) |>
  dplyr::mutate(
    category = factor(
      category,
      levels = c("1 arm", "2 arms", "3 arms", "4 arms", "5 arms")
    )
  )

bar_plot <- ggplot(
  summary_df,
  aes(x = category, y = pct, fill = color)
) +
  geom_col(
    lwd = 0.5,
    col = "black",
    fill = summary_df$color,
    width = 0.6
  ) +
  scale_fill_identity() +
  labs(
    x = "CTCF site presence",
    y = "Percentage (%)"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "gray70", linewidth = 0.3),
    plot.background = NULL,
    legend.position = "none",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )

ggsave(
  fs::path(plots, "ctcf_acro_arm_distribution.pdf"),
  bar_plot,
  width = 3.7,
  height = 1.25,
  bg = "transparent"
)

shared_sites <- odgi_pos_res |>
  dplyr::filter(
    !(
      query_chrom == ref_chrom &
        query_start == ref_start &
        query_end == ref_end
    )
  ) |>
  dplyr::left_join(
    stats_per_site,
    by = c("query_chrom", "query_start", "query_end")
  ) |>
  dplyr::mutate(
    color = adjustcolor(
      dplyr::recode(
        category,
        "1 arm"  = "#FCBBA1",
        "2 arms" = "#FB6A4A",
        "3 arms" = "#DE2D26",
        "4 arms" = "#A50F15",
        "5 arms" = "#67000D",
        .default = "grey90"
      ),
      alpha.f = 0.1
    )
  ) |>
  arrange(category)

shown_chrs <- sort(
  unique(c(shared_sites$query_chrom, shared_sites$ref_chrom))
)

censat_df <- as.data.frame(
  subsetByOverlaps(
    chm13_censat_gr,
    acrocentric_regions,
    ignore.strand = TRUE
  )
) |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::filter(chrom %in% shown_chrs) |>
  dplyr::mutate(chrom = as.character(chrom))

acro_xlim <- as.data.frame(acrocentric_regions) |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::mutate(chrom = as.character(chrom)) |>
  dplyr::group_by(chrom) |>
  dplyr::summarise(
    xmin = min(start),
    xmax = max(end),
    .groups = "drop"
  )

chm13_cytobands <- fs::path(data, "chm13", "chm13v2.0_cytobands_allchrs.bed")

chm13_cytobands_df <- read.table(
  chm13_cytobands,
  col.names = c("chrom", "start", "end", "name", "gieStain")
) |>
  dplyr::filter(chrom %in% shown_chrs) |>
  dplyr::left_join(acro_xlim, by = "chrom") |>
  dplyr::filter(start < xmax, end > xmin) |>
  dplyr::mutate(
    start = pmax(start, xmin),
    end   = pmin(end, xmax)
  ) |>
  dplyr::select(chrom, start, end, name, gieStain)

col_fun <- colorRamp2(
  c(
    0.9, 0.93, 0.98, 1
  ),
  c(
    "#66006F",
    "#6D6D6D",
    "#FCFF00",
    "#F29200"
  )
)

suppressMessages({
  chm13_sd_full <- read_tsv(
    fs::path(
      data, "chm13", "chm13v2.0_SD.full.bed"
    ),
    show_col_types = FALSE
  ) |>
    dplyr::select(
      `#chr1`, `start1...2`, `end1...3`,
      strand1, fracMatch
    ) |>
    dplyr::rename(
      seqnames = `#chr1`,
      start = `start1...2`,
      end = `end1...3`,
      strand = strand1,
      similarity = fracMatch
    )
})

chm13_sd_gr <- GRanges(
  seqnames = chm13_sd_full$seqnames,
  ranges   = IRanges(start = chm13_sd_full$start, end = chm13_sd_full$end),
  strand   = chm13_sd_full$strand,
  similarity = chm13_sd_full$similarity
)

acro_sds <- as.data.frame(
  subsetByOverlaps(
    chm13_sd_gr,
    acrocentric_regions,
    type = "within",
    ignore.strand = TRUE
  )
) |>
  dplyr::rename(chrom = seqnames) |>
  dplyr::mutate(itemRgb = col_fun(similarity))

chm13_acro_phrs <- as.data.frame(
  subsetByOverlaps(
    rtracklayer::import(
      fs::path(
        data, "chm13", "chm13_acro_PHRs.bed"
      )
    ),
    acrocentric_regions,
    type = "within",
    ignore.strand = TRUE
  )
) |>
  dplyr::rename(chrom = seqnames)

pdf(fs::path(plots, "acro_circos.pdf"), width = 4, height = 3.5)
lim <- 0.95

par(mar = c(0, 0, 0, 0))
circos.clear()
circos.par(
  "gap.degree" = 5,
  "start.degree" = 90,
  "cell.padding" = c(0, 0, 0, 0),
  "canvas.xlim" = c(-lim, lim),
  "canvas.ylim" = c(-lim, lim)
)

circos.initializeWithIdeogram(
  chm13_cytobands_df,
  chromosome.index = shown_chrs,
  plotType = c("axis", "labels")
)

circos.genomicTrack(
  dplyr::select(censat_df, c("chrom", "start", "end", "itemRgb")),
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region,
      ytop = 1,
      ybottom = 0,
      col  = value$itemRgb,
      border = NA
    )
  }
)

circos.genomicTrack(
  dplyr::select(acro_sds, c("chrom", "start", "end", "itemRgb")),
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region,
      ytop = 1,
      ybottom = 0,
      col = value$itemRgb,
      border = NA
    )
  }
)

circos.genomicTrack(
  dplyr::select(chm13_acro_phrs, c("chrom", "start", "end")),
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region,
      ytop = 1,
      ybottom = 0,
      col  = "#896949",
      border = NA
    )
  }
)

link_col <- shared_sites$color

circos.genomicLink(
  select(shared_sites, c(query_chrom, query_start, query_end)),
  select(shared_sites, c(ref_chrom, ref_start, ref_end)),
  col = link_col,
  border = link_col,
  lwd = 0.5
)

sd_legend <- Legend(
  at = c(
    0.9, 0.93, 0.98, 1
  ),
  col_fun = col_fun,
  title_position = "topleft",
  title = "SD Similarity",
  direction = "horizontal",
  labels_gp = gpar(fontsize = 5),
  title_gp  = gpar(fontsize = 6),
  legend_width  = unit(1.8, "cm"),
  grid_height   = unit(1.5, "mm")
)

draw(
  sd_legend,
  x = unit(0.95, "npc"),
  y = unit(0.04, "npc"),
  just = c("right", "bottom")
)

circos.clear()
dev.off()