#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript coverage_breadth_overlap.r <fivebase_sites_parquet> <ont_sites_parquet> <sample_id> <out_file>

How much of the genome each platform interrogates at the chosen min-
coverage threshold (from 2_persite_pileup/), and how much of that
overlaps between platforms - a basic breadth-of-coverage sanity check
alongside the depth-focused concordance metrics in 3_methylation_concordance/.
")
}

fivebase_sites_file <- args[[1]]
ont_sites_file      <- args[[2]]
sample_id           <- args[[3]]
out_file            <- args[[4]]

sites_5b  <- as.data.table(read_parquet(fivebase_sites_file))
sites_ont <- as.data.table(read_parquet(ont_sites_file))

n_shared <- nrow(merge(sites_5b[, .(chrom, cpg_pos)], sites_ont[, .(chrom, cpg_pos)], by = c("chrom", "cpg_pos")))

result <- data.table(
  sample_id = sample_id,
  n_sites_5base = nrow(sites_5b),
  n_sites_ont = nrow(sites_ont),
  n_sites_shared = n_shared,
  pct_5base_covered_by_ont = 100 * n_shared / nrow(sites_5b),
  pct_ont_covered_by_5base = 100 * n_shared / nrow(sites_ont)
)

fwrite_default(result, out_file)
print(result)
