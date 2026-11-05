# Load necessary libraries
suppressPackageStartupMessages({
  suppressWarnings({
    library("universalmotif", quietly = TRUE, verbose = FALSE)
    library("ggplot2", quietly = TRUE, verbose = FALSE)
  })
})

# Read the configuration YAML file that stores
# the paths of the directories used throughout the analysis
config <- read_yaml("../config.yml")

res_path <- fs::path(config$res) # Path containing all results

# Path containing metadata for T2T primates
data_path <- fs::path(config$data)
plots <- fs::path(res_path, "plots") # Path to store any generated plots

# Create a function that plots and saves each PWM
create_logo_plot <- function(
  input_file,
  plot_path
) {
  mot <- universalmotif::read_meme(input_file)

  p <- universalmotif::view_motifs(mot, use.type = "ICM") +
    theme(
      axis.text.x  = element_text(size = 4),
      axis.text.y  = element_text(size = 6),
      axis.title.y = element_text(size = 8),
      plot.background = NULL
    )

  ggplot2::ggsave(
    plot_path,
    p,
    width = 2.7,
    height = 1,
    bg = "transparent"
  )
}

# Detect input MEME files
meme_list <- list.files(
  path = fs::path(data_path, "motifs"),
  pattern = "\\.meme$",
  full.names = TRUE
)

# Iterate over input MEME files, create output plot paths, and create the
for (meme_file in meme_list) {
  out_name <- fs::path_ext_set(fs::path_file(meme_file), "pdf")
  create_logo_plot(
    input_file = meme_file,
    plot_path  = fs::path(plots, "figure_1", out_name)
  )
}