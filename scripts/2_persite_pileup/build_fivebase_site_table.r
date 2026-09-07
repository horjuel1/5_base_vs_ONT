#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript build_fivebase_site_table.r <mod_pileup_file> <out_parquet> [min_coverage]

Collapses 5-base's per-strand mod_pileup.txt (one row per (position,
strand)) into one row per CpG dinucleotide, so it can be joined against
the ONT per-site table on (chrom, cpg_pos).

mod_pileup.txt columns (tab-separated, first column literally '#chr'):
  #chr start end name beta_est strand unmod mod no_snp snp coverage
  genotype gt_p_score gt_conf_score
A CpG's '+' strand row (the C) sits at `start`; its '-' strand row (the
G, i.e. the C on the opposite strand) sits at `start+1` in this file's
convention, so shifting the '-' strand row back by 1bp aligns both
strands' counts onto the same cpg_pos before summing.

genotype != C/C (on +) or != G/G (on -) means the reference C is
actually a SNP in this sample - apparent (un)methylation there is a
genotyping artifact, not biology - so those positions are dropped by
default (require_no_snp).
")
}

pileup_file  <- args[[1]]
out_file     <- args[[2]]
min_coverage <- if (length(args) >= 3) as.integer(args[[3]]) else 5L
require_no_snp <- TRUE

if (!file.exists(pileup_file)) stop("mod_pileup_file not found: ", pileup_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", pileup_file)
dt <- fread(pileup_file, sep = "\t",
            select = c("#chr", "start", "strand", "unmod", "mod", "genotype"),
            col.names = c("chrom", "start", "strand", "unmod", "mod", "genotype"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt <- dt[chrom %in% CHROM_ORDER]

if (require_no_snp) {
  n_before <- nrow(dt)
  dt <- dt[(strand == "+" & genotype == "C/C") | (strand == "-" & genotype == "G/G")]
  log_msg("Positions after SNP/genotype filter: ", formatC(nrow(dt), format = "d", big.mark = ","),
          " / ", formatC(n_before, format = "d", big.mark = ","))
}

dt[, cpg_pos := ifelse(strand == "+", start, start - 1L)]

site_dt <- dt[, .(mod_count = sum(mod), coverage = sum(mod) + sum(unmod)), by = .(chrom, cpg_pos)]
site_dt <- site_dt[coverage >= min_coverage]
site_dt[, beta := mod_count / coverage]
log_msg("CpG sites with coverage >= ", min_coverage, ": ", formatC(nrow(site_dt), format = "d", big.mark = ","))

write_parquet(site_dt, out_file)
log_msg("Saved: ", out_file)
