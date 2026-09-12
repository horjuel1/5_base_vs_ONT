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

# Shared ONT calls.tsv filtering/classification, used by
# 2_persite_pileup/build_ont_site_table.r,
# 4_cpg_correlation_distance/2_compute_ont_cpg_pair_correlation.r, and
# 6_smmethid_compat/build_ont_moleculebed.r - previously duplicated across
# the first two verbatim, which is exactly how the missing secondary/
# supplementary-alignment filter would have silently diverged between them
# if fixed in only one place.
#
# `dt` must have columns: chrom, ref_position, ref_mod_strand, call_code,
# fail, within_alignment, flag (plus whatever else the caller wants kept,
# e.g. read_id, call_prob - passed through untouched).
#
# Filtering order matches smMethID's bin/calls_to_moleculebed.r: drop
# secondary/supplementary alignment records (SAM flag 0x100/0x800) first -
# a chimeric/concatemer ONT read's supplementary segment can call CpGs on a
# physically unrelated locus under the same read_id - then restrict to
# autosomes/within-alignment/passing calls, then classify call_code's first
# character as modified ('m'=5mC, 'h'=5hmC - TAPS/5-base's chemistry
# converts both, so they're pooled as one 'modified' state for a like-for-
# like cross-platform comparison) vs unmodified ('-'). Adds `is_mod` (0/1)
# and `cpg_pos` (both strands' CpG anchored to the same coordinate: '-'
# strand calls land on the CpG's G, shifted back 1bp to match the '+'
# strand C, same convention build_fivebase_site_table.r uses).
filter_and_classify_ont_calls <- function(dt, min_call_prob = NA_real_) {
  dt <- data.table::copy(dt)
  dt[, flag := as.integer(flag)]
  n_before_flag <- nrow(dt)
  dt <- dt[bitwAnd(flag, 0x900L) == 0L]
  log_msg("Rows from primary alignments (flag & 0x900 == 0): ", formatC(nrow(dt), format = "d", big.mark = ","),
          " / ", formatC(n_before_flag, format = "d", big.mark = ","))

  dt <- dt[chrom %in% CHROM_ORDER & within_alignment == TRUE & fail == FALSE]
  mod_char <- substr(dt$call_code, 1, 1)
  dt <- dt[mod_char %in% c("m", "h", "-")]
  mod_char <- mod_char[mod_char %in% c("m", "h", "-")]
  if (!is.na(min_call_prob)) {
    keep <- dt$call_prob >= min_call_prob
    dt <- dt[keep]
    mod_char <- mod_char[keep]
  }
  log_msg("Confident calls after filtering: ", formatC(nrow(dt), format = "d", big.mark = ","))

  dt[, is_mod := as.integer(mod_char %in% c("m", "h"))]
  dt[, cpg_pos := ifelse(ref_mod_strand == "-", ref_position - 1L, ref_position)]
  dt[]
}
