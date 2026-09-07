#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("
Usage:
Rscript methylation_by_fragment_length.r <fivebase_mod_frag_file> <ont_cpg_density_file> <sample_id> <out_tsv> <out_png>

Does apparent methylation / CpG density depend on fragment length -
relevant for cfDNA, where short fragments are enriched at nucleosome-
protected, often CpG-dense regions (a known DELFI/fragmentomics-
adjacent signal, not just a QC nicety here).

ont_cpg_density_file columns (tab-separated, gzipped):
  read_id chrom fragment_length genomic_cpg_count called_cpg_count methylated_cpg_count
")
}

fivebase_frag_file    <- args[[1]]
ont_cpg_density_file  <- args[[2]]
sample_id             <- args[[3]]
out_tsv               <- args[[4]]
out_png               <- args[[5]]

LENGTH_BREAKS <- c(0, 100, 150, 167, 200, 250, 320, 400, 600, 1000, 5000, Inf)

log_msg("Reading 5-base: ", fivebase_frag_file)
dt_5b <- fread(fivebase_frag_file, sep = "\t",
               select = c("#chr", "read_id", "insert_size", "num_cpg", "num_mod"),
               col.names = c("chrom", "read_id", "insert_size", "num_cpg", "num_mod"))
dt_5b <- dt_5b[, .(fragment_length = insert_size[1], num_cpg = sum(num_cpg), num_mod = sum(num_mod)), by = read_id]
dt_5b <- dt_5b[num_cpg > 0]
dt_5b[, meth_frac := num_mod / num_cpg]
dt_5b[, cpg_density := num_cpg / fragment_length]

log_msg("Reading ONT: ", ont_cpg_density_file)
dt_ont <- fread(ont_cpg_density_file, sep = "\t")
dt_ont <- dt_ont[called_cpg_count > 0]
dt_ont[, meth_frac := methylated_cpg_count / called_cpg_count]
dt_ont[, cpg_density := genomic_cpg_count / fragment_length]

bin_summary <- function(dt, label) {
  dt <- copy(dt)
  # bin_idx (an integer into LENGTH_BREAKS) instead of cut()'s factor label,
  # so the plotting midpoint below comes from the actual break values, not
  # from parsing strings back out of a label like "[5e+03,Inf)" - a regex
  # over that label would silently mis-parse scientific notation.
  dt[, bin_idx := findInterval(fragment_length, LENGTH_BREAKS, rightmost.closed = FALSE)]
  g <- dt[, .(n = .N, mean_meth_frac = mean(meth_frac), mean_cpg_density = mean(cpg_density)), by = bin_idx]
  g[, platform := label]
  g
}

summary_dt <- rbind(bin_summary(dt_5b, "5-base"), bin_summary(dt_ont, "ONT"))
summary_dt[, length_bin := sprintf("[%s,%s)", LENGTH_BREAKS[bin_idx], LENGTH_BREAKS[bin_idx + 1])]
summary_dt[, sample_id := sample_id]
fwrite_default(summary_dt[, .(sample_id, platform, length_bin, n, mean_meth_frac, mean_cpg_density)], out_tsv)
print(summary_dt)

# Midpoint of each bin for plotting on a continuous x-axis; the final
# open-ended bin ([5000,Inf)) has no finite midpoint, so it's dropped from
# the plot (still present in the TSV).
summary_dt[, x := (LENGTH_BREAKS[bin_idx] + LENGTH_BREAKS[bin_idx + 1]) / 2]
plot_dt <- summary_dt[is.finite(x)]

p1 <- ggplot(plot_dt, aes(x = x, y = mean_meth_frac, color = platform)) +
  geom_line(linewidth = 0.8) + geom_point() +
  scale_color_manual(values = c("5-base" = PALETTE$fivebase, "ONT" = PALETTE$ont)) +
  labs(x = "Fragment length (bp)", y = "Mean methylation fraction", title = "Methylation vs fragment length") +
  plot_theme()

p2 <- ggplot(plot_dt, aes(x = x, y = mean_cpg_density, color = platform)) +
  geom_line(linewidth = 0.8) + geom_point() +
  scale_color_manual(values = c("5-base" = PALETTE$fivebase, "ONT" = PALETTE$ont)) +
  labs(x = "Fragment length (bp)", y = "CpGs per bp", title = "CpG density vs fragment length") +
  plot_theme()

ggsave_default(out_png, p1 + p2, width = 12, height = 5)
log_msg("Done.")
