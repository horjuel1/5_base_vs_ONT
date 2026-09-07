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
Rscript build_ont_site_table.r <ont_calls_file> <out_parquet> [min_coverage] [min_call_prob]

Aggregates ONT's per-call table (one row per CpG call per read) into one
row per CpG site, analogous to a modkit-pileup bedMethyl, so it can be
joined against the 5-base per-site table on (chrom, cpg_pos).

calls.tsv columns used here (tab-separated, header present):
  read_id ref_position chrom ref_mod_strand call_prob call_code
  fail within_alignment
call_code is assumed 'm'=5mC, 'h'=5hmC, '-'=canonical/unmodified (first
character). 5-base's TAPS chemistry converts both 5mC and 5hmC, so m+h
are counted together as 'modified' here for a like-for-like comparison.
`fail` is assumed to already flag low-confidence/ambiguous calls to
exclude - confirm this against the calling pipeline that produced
calls.tsv before trusting downstream numbers; min_call_prob is available
as an extra, independent filter on top of `fail` if that assumption
turns out to be wrong or too permissive.

ref_position is assumed to be the position of the called base on the
reference in the direction given by ref_mod_strand; for '-' strand calls
that position is the CpG's G, so it's shifted back by 1bp to land on the
same cpg_pos convention as the '+' strand C (matching
build_fivebase_site_table.r's convention).
")
}

calls_file    <- args[[1]]
out_file      <- args[[2]]
min_coverage  <- if (length(args) >= 3) as.integer(args[[3]]) else 5L
min_call_prob <- if (length(args) >= 4) as.numeric(args[[4]]) else NA_real_

if (!file.exists(calls_file)) stop("ont_calls_file not found: ", calls_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", calls_file)
dt <- fread(calls_file, sep = "\t",
            select = c("chrom", "ref_position", "ref_mod_strand", "call_prob",
                       "call_code", "fail", "within_alignment"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt <- dt[chrom %in% CHROM_ORDER & within_alignment == TRUE & fail == FALSE]
mod_char <- substr(dt$call_code, 1, 1)
dt <- dt[mod_char %in% c("m", "h", "-")]
mod_char <- mod_char[mod_char %in% c("m", "h", "-")]
if (!is.na(min_call_prob)) dt <- dt[call_prob >= min_call_prob]
log_msg("Confident calls after filtering: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt[, is_mod := as.integer(mod_char %in% c("m", "h"))]
dt[, cpg_pos := ifelse(ref_mod_strand == "-", ref_position - 1L, ref_position)]

site_dt <- dt[, .(mod_count = sum(is_mod), coverage = .N), by = .(chrom, cpg_pos)]
site_dt <- site_dt[coverage >= min_coverage]
site_dt[, beta := mod_count / coverage]
log_msg("CpG sites with coverage >= ", min_coverage, ": ", formatC(nrow(site_dt), format = "d", big.mark = ","))

write_parquet(site_dt, out_file)
log_msg("Saved: ", out_file)
