#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript 3_plot_correlation_by_distance_comparison.r <fivebase_result_file> <ont_result_file> <sample_id> <out_png>

Overlays the two platforms' CpG-pair correlation-by-distance curves
(outputs of 1_compute_fivebase_... and 2_compute_ont_...).
")
}

fivebase_file <- args[[1]]
ont_file      <- args[[2]]
sample_id     <- args[[3]]
out_png       <- args[[4]]

dt <- rbind(fread(fivebase_file), fread(ont_file))

p <- ggplot(dt, aes(x = bin_mid, y = mean_phi, color = platform)) +
  geom_line(aes(group = platform), linewidth = 0.6) +
  geom_point(aes(size = n_pairs), alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("5-base" = PALETTE$fivebase, "ONT" = PALETTE$ont)) +
  scale_size_continuous(name = "N CpG pairs", trans = "log10") +
  labs(x = "Genomic distance between CpG pair (bp, log scale)",
       y = "Within-read/fragment methylation correlation (mean phi)",
       title = paste0(sample_id, ": CpG co-methylation correlation vs. distance")) +
  plot_theme()

ggsave_default(out_png, p, width = 8, height = 5.5)
log_msg("Done.")
