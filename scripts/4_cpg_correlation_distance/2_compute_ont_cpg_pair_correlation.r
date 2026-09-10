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
Rscript 2_compute_ont_cpg_pair_correlation.r <ont_calls_file> <sample_id> <out_file> [max_distance_bp]

Within-read CpG co-methylation correlation as a function of genomic
distance, for ONT. calls.tsv is already one row per (read_id, CpG call),
so unlike the 5-base side there's no comma-list to explode - it's
already in the long (read_id, chrom, pos, state) format the shared
engine (_cpg_pair_correlation_utils.r) expects.

Reads blacklisted regions via PlasmaTools::filters.hg19 if that package
is installed (cluster-only; skipped with a warning if not available -
e.g. when testing locally) - repeat-mapping ONT reads piling into
satellite DNA/centromeres/rDNA are a documented OOM risk for this exact
bucketed-correlation approach on the buffy coat side of this lab's other
project, so filtering them out here too is cheap insurance, not
speculative.

Pass a much larger max_distance_bp than the 5-base run to see how far
long reads let this correlation persist - short fragments/reads
fundamentally can't be measured past their own length, so this is where
ONT's long-read advantage should show up most clearly.
")
}

calls_file       <- args[[1]]
sample_id        <- args[[2]]
out_file         <- args[[3]]
max_distance_bp  <- if (length(args) >= 4) as.integer(args[[4]]) else 2000L

if (!file.exists(calls_file)) stop("ont_calls_file not found: ", calls_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

BIN_WIDTH_BP <- 10L
BUCKET_SIZE_BP <- 20000L
MIN_READS_PER_PAIR <- 5L
MAX_READS_PER_BUCKET <- 20000L
MAX_CPGS_PER_BUCKET  <- 3000L

if (max_distance_bp >= BUCKET_SIZE_BP) {
  stop("max_distance_bp (", max_distance_bp, ") must be < BUCKET_SIZE_BP (", BUCKET_SIZE_BP,
       ") for the core+halo bucketing to be correct - raise BUCKET_SIZE_BP in this script if you need a wider distance range.")
}

log_msg("Reading: ", calls_file)
dt <- fread(calls_file, sep = "\t",
            select = c("read_id", "chrom", "ref_position", "ref_mod_strand", "call_code", "fail", "within_alignment", "flag"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

# Filtering (including dropping secondary/supplementary alignment calls -
# see filter_and_classify_ont_calls()'s own comment for why that matters
# more here than in 2_persite_pileup/build_ont_site_table.r: a chimeric/
# concatemer read's supplementary segment would otherwise fabricate within-
# read CpG pairs across two physically unrelated loci) is shared with
# build_ont_site_table.r via _shared_utils.r, so both stay in sync.
dt <- filter_and_classify_ont_calls(dt)
dt[, state := is_mod]
dt[, pos := cpg_pos]
long_dt <- unique(dt[, .(read_id, chrom, pos, state)], by = c("read_id", "chrom", "pos"))
log_msg("Read-level CpG calls: ", formatC(nrow(long_dt), format = "d", big.mark = ","))

blacklist_dt <- NULL
if (requireNamespace("PlasmaTools", quietly = TRUE)) {
  log_msg("Loading blacklist (PlasmaTools::filters.hg19)...")
  bl <- as.data.table(PlasmaTools::filters.hg19)
  if ("seqnames" %in% names(bl)) setnames(bl, "seqnames", "chrom")
  blacklist_dt <- bl[, .(chrom = as.character(chrom), start = as.integer(start), end = as.integer(end))]
  blacklist_dt <- blacklist_dt[chrom %in% CHROM_ORDER & !is.na(start) & !is.na(end) & start < end]
  blacklist_dt <- unique(blacklist_dt)
} else {
  log_msg("PlasmaTools not installed - skipping blacklist filtering (fine for local testing; ",
          "install PlasmaTools on the cluster for a production run).")
}

result <- compute_pair_correlation_by_distance(
  long_dt, max_distance_bp = max_distance_bp, bin_width_bp = BIN_WIDTH_BP,
  bucket_size_bp = BUCKET_SIZE_BP, min_reads_per_pair = MIN_READS_PER_PAIR,
  max_reads_per_bucket = MAX_READS_PER_BUCKET, max_cpgs_per_bucket = MAX_CPGS_PER_BUCKET,
  blacklist_dt = blacklist_dt
)
result[, platform := "ONT"]
result[, sample_id := sample_id]

fwrite_default(result, out_file)
print(result)
log_msg("\nDone.")
