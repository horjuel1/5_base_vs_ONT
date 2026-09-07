#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "_peak_detection_utils.r"))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript compute_fragment_length_distribution.r <fivebase_mod_frag_file> <ont_fraglen_file> <sample_id> <out_dir>

Fragment length distribution, mono-/di-nucleosome periodicity, and
longer-fragment abundance for one sample sequenced on both platforms.

fivebase_mod_frag_file (5-base/TAPS, Illumina short-read; tab-separated,
first column literally '#chr'):
  #chr start end read_id mapq orientation insert_size read_length flag
  num_cpg num_mod mod_cpgs unmod_cpgs snp_cpgs
Each PHYSICAL fragment appears as TWO rows (one per mate), same read_id
and insert_size, opposite orientation - must dedupe by read_id or every
fragment is double-counted (same gotcha as
buffy_methylation_architecture/scripts/cfdna_nucleosome_periodicity/
1_compute_fragment_length_distribution.r, which this script's 5-base
histogram/peak logic is adapted from).

ont_fraglen_file (ONT long-read; tab-separated, gzipped):
  read_id chrom fragment_length read_length
Already one row per read - no dedup needed.

Unlike the short-read-only cfDNA periodicity script this is adapted
from, ONT fragments routinely exceed the 50-1000bp short-read-only
filtering window (long native DNA reads, no size-selection/PCR
artifact). So ONT gets its own unclipped summary + a wide-range
long-fragment abundance table, while BOTH platforms additionally get a
shared 50-1000bp/10bp-bin histogram (mono-/di-nucleosome periodicity
range) so the two platforms are visually and numerically comparable on
that shared axis.
")
}

fivebase_frag_file <- args[[1]]
ont_fraglen_file    <- args[[2]]
sample_id           <- args[[3]]
out_dir             <- args[[4]]

if (!file.exists(fivebase_frag_file)) stop("fivebase_mod_frag_file not found: ", fivebase_frag_file)
if (!file.exists(ont_fraglen_file))   stop("ont_fraglen_file not found: ", ont_fraglen_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────
CORE_MIN_BP   <- 50    # below this is almost always adapter dimer/artifact, not a real fragment
CORE_MAX_BP   <- 1000  # mono- through ~5-nucleosome-repeat range shared by both platforms
BIN_WIDTH_BP  <- 10    # matches buffy_methylation_architecture's periodicity bin width
SMOOTH_WINDOW_BINS <- 3
MIN_PEAK_PROMINENCE_FRAC <- 0.02
MIN_SPACING_PEAK_PROMINENCE_FRAC <- 0.05
N_SPACING_PEAKS <- 5

# Longer-fragment abundance cutoffs: mono-nucleosome, standard cfDNA
# "long", di-nucleosome, high-molecular-weight/leukocyte-contamination scale.
LONG_CUTOFFS_BP <- c(167, 250, 500, 1000, 5000, 10000)

# DELFI-style short:long ratio (Cristiano et al.) - a standard cfDNA
# fragmentomics QC metric, directly relevant since this is a LUCAS cfDNA sample.
SHORT_RANGE <- c(100, 150)
LONG_RANGE  <- c(151, 220)

# ─────────────────────────────────────────────
# Load
# ─────────────────────────────────────────────
log_msg("Reading 5-base: ", fivebase_frag_file)
dt_5b <- fread(fivebase_frag_file, sep = "\t", select = c("#chr", "read_id", "insert_size"),
               col.names = c("chr", "read_id", "insert_size"))
dt_5b <- unique(dt_5b[, .(read_id, fragment_length = insert_size)])
dt_5b <- dt_5b[!is.na(fragment_length) & fragment_length > 0]
log_msg("5-base unique fragments: ", formatC(nrow(dt_5b), format = "d", big.mark = ","))

log_msg("Reading ONT: ", ont_fraglen_file)
dt_ont <- fread(ont_fraglen_file, sep = "\t", select = c("read_id", "fragment_length"))
dt_ont <- dt_ont[!is.na(fragment_length) & fragment_length > 0]
log_msg("ONT reads: ", formatC(nrow(dt_ont), format = "d", big.mark = ","))

lengths_list <- list("5-base" = dt_5b$fragment_length, "ONT" = dt_ont$fragment_length)

# ─────────────────────────────────────────────
# Full-range summary + longer-fragment abundance (no upper clip - this is
# the whole point of comparing a platform whose reads can be tens of kb
# against one capped at ~short-read insert sizes).
# ─────────────────────────────────────────────
summarize_lengths <- function(x) {
  pct <- quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
  data.table(
    n = length(x), mean = mean(x), sd = sd(x),
    p5 = pct[[1]], p25 = pct[[2]], median = pct[[3]], p75 = pct[[4]], p95 = pct[[5]],
    max = max(x)
  )
}

long_abundance <- function(x) {
  as.data.table(as.list(setNames(
    sapply(LONG_CUTOFFS_BP, function(c) mean(x > c) * 100),
    paste0("pct_gt_", LONG_CUTOFFS_BP, "bp")
  )))
}

short_long_ratio <- function(x) {
  short_n <- sum(x >= SHORT_RANGE[1] & x <= SHORT_RANGE[2])
  long_n  <- sum(x >= LONG_RANGE[1] & x <= LONG_RANGE[2])
  if (long_n == 0) return(NA_real_)
  short_n / long_n
}

summary_dt <- rbindlist(lapply(names(lengths_list), function(label) {
  x <- lengths_list[[label]]
  cbind(data.table(platform = label), summarize_lengths(x), long_abundance(x),
        data.table(short_long_ratio = short_long_ratio(x)))
}))
summary_dt[, sample_id := sample_id]
fwrite_default(summary_dt, file.path(out_dir, "fragment_length_summary.tsv"))
print(summary_dt)

# ─────────────────────────────────────────────
# Shared core-range (50-1000bp) histogram + periodicity peak detection,
# one platform at a time, using buffy_methylation_architecture's already-
# verified smoothing/peak-finding method identically on both sides (see
# _peak_detection_utils.r header for why symmetry between sides matters).
# ─────────────────────────────────────────────
build_core_histogram <- function(x, label) {
  x <- x[x >= CORE_MIN_BP & x <= CORE_MAX_BP]
  dt <- data.table(length_bin = (x %/% BIN_WIDTH_BP) * BIN_WIDTH_BP)
  hist_dt <- dt[, .(n_fragments = .N), by = length_bin]
  setorder(hist_dt, length_bin)
  full_bins <- data.table(length_bin = seq(CORE_MIN_BP, CORE_MAX_BP, by = BIN_WIDTH_BP))
  hist_dt <- merge(full_bins, hist_dt, by = "length_bin", all.x = TRUE)
  hist_dt[is.na(n_fragments), n_fragments := 0L]
  hist_dt[, density := n_fragments / sum(n_fragments)]
  hist_dt[, density_smooth := smooth_signal(density, SMOOTH_WINDOW_BINS)]

  min_prominence <- MIN_PEAK_PROMINENCE_FRAC * max(hist_dt$density_smooth)
  hist_dt[, is_peak := find_local_maxima(density_smooth, min_prominence)]

  peaks_dt <- hist_dt[is_peak == TRUE][order(-density_smooth)]
  min_spacing_prominence <- MIN_SPACING_PEAK_PROMINENCE_FRAC * max(hist_dt$density_smooth)
  spacing_result <- spacing_from_top_peaks(peaks_dt, "length_bin", "density_smooth",
                                            min_spacing_prominence = min_spacing_prominence,
                                            n_spacing_peaks = N_SPACING_PEAKS)
  if (length(spacing_result$spacings) > 0) {
    log_msg(label, " peak spacing (bp): ", paste(spacing_result$spacings, collapse = ", "),
            " | mean: ", round(spacing_result$mean_spacing, 1))
  } else {
    log_msg(label, ": not enough eligible peaks to estimate spacing")
  }

  hist_dt[, platform := label]
  peaks_dt[, platform := label]
  list(hist = hist_dt, peaks = peaks_dt, mean_spacing = spacing_result$mean_spacing)
}

core_5b  <- build_core_histogram(lengths_list[["5-base"]], "5-base")
core_ont <- build_core_histogram(lengths_list[["ONT"]], "ONT")

hist_combined  <- rbind(core_5b$hist, core_ont$hist)
peaks_combined <- rbind(core_5b$peaks, core_ont$peaks)
hist_combined[, sample_id := sample_id]
peaks_combined[, sample_id := sample_id]

fwrite_default(hist_combined, file.path(out_dir, "fragment_length_core_histogram.tsv"))
fwrite_default(peaks_combined, file.path(out_dir, "fragment_length_core_peaks.tsv"))

# ─────────────────────────────────────────────
# Plots: core-range histogram overlay + full-range log-scale ECDF.
# ─────────────────────────────────────────────
p_hist <- ggplot(hist_combined, aes(x = length_bin, y = density_smooth, color = platform)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = peaks_combined, aes(x = length_bin, y = density_smooth), size = 2) +
  scale_color_manual(values = c("5-base" = PALETTE$fivebase, "ONT" = PALETTE$ont)) +
  labs(x = "Fragment length (bp)", y = "Smoothed density",
       title = paste0(sample_id, ": fragment length distribution (", CORE_MIN_BP, "-", CORE_MAX_BP, "bp)")) +
  plot_theme()

ecdf_dt <- rbindlist(lapply(names(lengths_list), function(label) {
  x <- sort(lengths_list[[label]])
  data.table(platform = label, length = x, ecdf = seq_len(length(x)) / length(x))
}))
p_ecdf <- ggplot(ecdf_dt, aes(x = length, y = ecdf, color = platform)) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  scale_color_manual(values = c("5-base" = PALETTE$fivebase, "ONT" = PALETTE$ont)) +
  labs(x = "Fragment length (bp, log scale)", y = "Cumulative fraction",
       title = paste0(sample_id, ": fragment length ECDF (full range)")) +
  plot_theme()

combined_plot <- patchwork::wrap_plots(p_hist, p_ecdf, nrow = 1)
ggsave_default(file.path(out_dir, "fragment_length_distribution.png"), combined_plot, width = 12, height = 5)

log_msg("\nDone.")
