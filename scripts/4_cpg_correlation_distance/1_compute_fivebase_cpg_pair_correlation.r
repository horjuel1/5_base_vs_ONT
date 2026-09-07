#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))
source(file.path(script_dir, "_cpg_pair_correlation_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript 1_compute_fivebase_cpg_pair_correlation.r <mod_frag_file> <sample_id> <out_file> [max_distance_bp]

Within-fragment CpG co-methylation correlation as a function of genomic
distance, for 5-base. mod_cpgs/unmod_cpgs are comma-separated 0-based
offsets from `start`; both mates of a read_id are pooled into one
fragment-level CpG list before pairing, since together they reconstruct
the original fragment's full methylation pattern (a short-read
fragment's own reach is inherently capped near its insert size - this is
expected to look very different from ONT's curve, not a bug in either).
See _cpg_pair_correlation_utils.r for the correlation/bucketing engine,
shared with the ONT side so both curves are computed identically.
")
}

frag_file        <- args[[1]]
sample_id        <- args[[2]]
out_file         <- args[[3]]
max_distance_bp  <- if (length(args) >= 4) as.integer(args[[4]]) else 2000L

if (!file.exists(frag_file)) stop("mod_frag_file not found: ", frag_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

BIN_WIDTH_BP <- 10L
BUCKET_SIZE_BP <- 20000L
MIN_READS_PER_PAIR <- 5L        # thinner cfDNA coverage than the buffy-coat side's MIN_READS_PER_PAIR=10
MAX_READS_PER_BUCKET <- 20000L
MAX_CPGS_PER_BUCKET  <- 3000L

log_msg("Reading: ", frag_file)
dt <- fread(frag_file, sep = "\t", select = c("#chr", "start", "read_id", "mod_cpgs", "unmod_cpgs"),
            col.names = c("chrom", "start", "read_id", "mod_cpgs", "unmod_cpgs"))
dt <- dt[chrom %in% CHROM_ORDER]
log_msg("Rows: ", formatC(nrow(dt), format = "d", big.mark = ","))

# Explode comma-separated offsets into one row per (read_id, chrom, pos,
# state). strsplit("", ",") returns character(0), so reads with zero
# mod (or zero unmod) CpGs correctly contribute zero rows here - no NA
# offsets to filter out, unlike a naive split-on-empty-string approach.
build_offset_rows <- function(dt, col, state_value) {
  split_list <- strsplit(dt[[col]], ",", fixed = TRUE)
  n <- lengths(split_list)
  offsets <- suppressWarnings(as.integer(unlist(split_list)))
  data.table(
    read_id = rep(dt$read_id, n),
    chrom = rep(dt$chrom, n),
    pos = rep(dt$start, n) + offsets,
    state = state_value
  )
}

long_dt <- rbind(build_offset_rows(dt, "mod_cpgs", 1L), build_offset_rows(dt, "unmod_cpgs", 0L))
long_dt <- long_dt[!is.na(pos)]
# Mates of the same fragment can both cover the same CpG in an overlapping
# insert - collapse to one row per (read_id, pos) before pairing so that
# position isn't counted twice within a single fragment's CpG list.
long_dt <- unique(long_dt, by = c("read_id", "chrom", "pos"))
log_msg("Fragment-level CpG calls: ", formatC(nrow(long_dt), format = "d", big.mark = ","))

result <- compute_pair_correlation_by_distance(
  long_dt, max_distance_bp = max_distance_bp, bin_width_bp = BIN_WIDTH_BP,
  bucket_size_bp = BUCKET_SIZE_BP, min_reads_per_pair = MIN_READS_PER_PAIR,
  max_reads_per_bucket = MAX_READS_PER_BUCKET, max_cpgs_per_bucket = MAX_CPGS_PER_BUCKET
)
result[, platform := "5-base"]
result[, sample_id := sample_id]

fwrite_default(result, out_file)
print(result)
log_msg("\nDone.")
