#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(patchwork)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript compute_methylation_concordance.r <fivebase_sites_parquet> <ont_sites_parquet> <sample_id> <out_dir> [meth_threshold]

Joins the two platforms' per-CpG-site tables (from
2_persite_pileup/build_*_site_table.r) on (chrom, cpg_pos) and reports
per-CpG methylation concordance: Pearson/Spearman correlation, RMSD/MAD,
Bland-Altman bias + limits of agreement, binary agreement/sensitivity/
specificity/Cohen's kappa at a beta cutoff, and how all of the above
change with coverage depth.
")
}

fivebase_sites_file <- args[[1]]
ont_sites_file      <- args[[2]]
sample_id           <- args[[3]]
out_dir             <- args[[4]]
meth_threshold      <- if (length(args) >= 5) as.numeric(args[[5]]) else 0.5

if (!file.exists(fivebase_sites_file)) stop("fivebase_sites_parquet not found: ", fivebase_sites_file)
if (!file.exists(ont_sites_file))      stop("ont_sites_parquet not found: ", ont_sites_file)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading site tables...")
sites_5b  <- as.data.table(read_parquet(fivebase_sites_file))
sites_ont <- as.data.table(read_parquet(ont_sites_file))
setnames(sites_5b, c("mod_count", "coverage", "beta"), c("mod_count_5b", "cov_5base", "beta_5base"))
setnames(sites_ont, c("mod_count", "coverage", "beta"), c("mod_count_ont", "cov_ont", "beta_ont"))

joined <- merge(sites_5b, sites_ont, by = c("chrom", "cpg_pos"))
log_msg("Shared CpG sites: ", formatC(nrow(joined), format = "d", big.mark = ","))
if (nrow(joined) < 10) stop("Too few shared CpG sites to compute concordance (", nrow(joined), ") - check chrom naming/build match between platforms.")

write_parquet(joined, file.path(out_dir, "joined_cpg_sites.parquet"))
log_msg("Saved: ", file.path(out_dir, "joined_cpg_sites.parquet"))

# ─────────────────────────────────────────────
# Concordance metrics - one function, reused for the overall summary and
# for each coverage-stratified subset below.
# ─────────────────────────────────────────────
concordance_metrics <- function(dt, threshold) {
  x <- dt$beta_5base
  y <- dt$beta_ont
  diff <- x - y

  call_x <- x >= threshold
  call_y <- y >= threshold
  agree <- mean(call_x == call_y)
  tp <- sum(call_x & call_y); tn <- sum(!call_x & !call_y)
  fp <- sum(call_x & !call_y); fn <- sum(!call_x & call_y)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  po <- agree
  pe <- mean(call_x) * mean(call_y) + mean(!call_x) * mean(!call_y)
  kappa <- if (pe != 1) (po - pe) / (1 - pe) else NA_real_

  data.table(
    n_sites = nrow(dt),
    pearson_r = cor(x, y, method = "pearson"),
    spearman_r = cor(x, y, method = "spearman"),
    rmsd = sqrt(mean(diff^2)),
    mad = mean(abs(diff)),
    bias_mean_diff_5base_minus_ont = mean(diff),
    loa_lower = mean(diff) - 1.96 * sd(diff),
    loa_upper = mean(diff) + 1.96 * sd(diff),
    binary_agreement = agree,
    sensitivity_5base_vs_ont = sens,
    specificity_5base_vs_ont = spec,
    cohens_kappa = kappa,
    global_mean_beta_5base = weighted.mean(x, dt$cov_5base),
    global_mean_beta_ont = weighted.mean(y, dt$cov_ont)
  )
}

overall <- concordance_metrics(joined, meth_threshold)
overall[, sample_id := sample_id]
fwrite_default(overall, file.path(out_dir, "methylation_concordance_summary.tsv"))
print(overall)

# ─────────────────────────────────────────────
# Coverage-stratified concordance - does agreement improve with depth on
# both platforms.
# ─────────────────────────────────────────────
cov_edges <- c(0, 5, 10, 20, 50, Inf)
cov_labels <- c("0-5", "5-10", "10-20", "20-50", "50+")
joined[, min_cov := pmin(cov_5base, cov_ont)]
joined[, cov_bin := cut(min_cov, breaks = cov_edges, labels = cov_labels, right = FALSE)]

by_cov <- joined[, if (.N >= 10) concordance_metrics(.SD, meth_threshold), by = cov_bin]
by_cov[, sample_id := sample_id]
fwrite_default(by_cov, file.path(out_dir, "methylation_concordance_by_coverage.tsv"))
print(by_cov)

# ─────────────────────────────────────────────
# Plots: hexbin scatter + Bland-Altman.
# ─────────────────────────────────────────────
joined[, mean_beta := (beta_5base + beta_ont) / 2]
joined[, diff_beta := beta_5base - beta_ont]

p_scatter <- ggplot(joined, aes(x = beta_5base, y = beta_ont)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, color = PALETTE$reference, linetype = "dashed") +
  scale_fill_viridis_c(trans = "log10", name = "N sites") +
  labs(x = "5-base beta", y = "ONT beta",
       title = paste0(sample_id, ": per-CpG methylation concordance"),
       subtitle = sprintf("Pearson r = %.3f, n = %s sites", overall$pearson_r, formatC(overall$n_sites, format = "d", big.mark = ","))) +
  plot_theme()

p_ba <- ggplot(joined, aes(x = mean_beta, y = diff_beta)) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_hline(yintercept = overall$bias_mean_diff_5base_minus_ont, color = "red", linetype = "dashed") +
  geom_hline(yintercept = overall$loa_lower, color = PALETTE$reference, linetype = "dotted") +
  geom_hline(yintercept = overall$loa_upper, color = PALETTE$reference, linetype = "dotted") +
  labs(x = "Mean beta (5-base, ONT)", y = "Beta difference (5-base - ONT)",
       title = "Bland-Altman") +
  plot_theme()

combined_plot <- p_scatter + p_ba
ggsave_default(file.path(out_dir, "methylation_concordance.png"), combined_plot, width = 12, height = 5)

log_msg("\nDone.")
