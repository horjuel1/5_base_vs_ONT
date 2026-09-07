# Sourced (not run standalone) by scripts across every topic folder in this
# project. Mirrors the sourced-utils convention used in
# buffy_methylation_architecture/scripts/cfdna_nucleosome_periodicity/
# (_peak_detection_utils.r) - one copy of logging/plotting/constants shared
# by both platforms' scripts, so neither side drifts from the other.

log_msg <- function(...) {
  message(...)
  flush(stderr())
}

# chr1-22 only (autosomes) - avoids sex-chromosome coverage/dosage
# differences confounding cross-platform methylation comparisons. Matches
# buffy_methylation_architecture's own chrom_order convention (though that
# project also keeps chrX; we drop it here since it's not load-bearing for
# a single-sample fragmentomics/methylation comparison).
CHROM_ORDER <- paste0("chr", 1:22)

# Hex palette matching buffy_methylation_architecture's plots, so figures
# from this project read as part of the same body of work.
PALETTE <- list(
  fivebase   = "#2a78d6",  # blue
  ont        = "#e34948",  # red
  fivebase_dark = "#184f95",
  ont_dark   = "#a02322",
  reference  = "#898781"   # gray, for y=x / threshold lines
)

plot_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size)
}

ggsave_default <- function(path, plot, width = 7, height = 5, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  log_msg("Saved: ", path)
}

fwrite_default <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, path, sep = "\t")
  log_msg("Saved: ", path)
}
