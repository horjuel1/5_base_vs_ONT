# Sourced (not run standalone) by 1_compute_fivebase_cpg_pair_correlation.r
# and 2_compute_ont_cpg_pair_correlation.r.
#
# Adapted from buffy_methylation_architecture/scripts/4_compute_region_metrics/
# 16_compute_cpg_pair_correlation_unconstrained_per_sample.r's core+halo
# bucketing/correlation engine, generalized to take an already-long
# (read_id, chrom, pos, state) table as input instead of parsing per-read
# comma-separated cpg_ids/calls strings - both platforms' raw files reduce
# to that same long format (see the two caller scripts), so one engine can
# serve both instead of duplicating/risking drift between two versions.
#
# Core + halo bucketing (unchanged rationale from the original): each
# bucket i owns a CORE range [i*BUCKET_SIZE_BP, (i+1)*BUCKET_SIZE_BP). A
# CpG at position p is gathered into bucket i's computation if
# i*BUCKET_SIZE_BP <= p < (i+1)*BUCKET_SIZE_BP + MAX_DISTANCE_BP - wide
# enough that any pair whose LOWER (left) CpG falls in this bucket's core
# is fully supported, even if its RIGHT CpG extends past the core into the
# next bucket. After computing every pair in the gathered set, only pairs
# whose lower position falls inside this bucket's own core are kept - so
# each physical pair is attributed to exactly one bucket, never double-
# counted, and never missing its full supporting read set.
#
# Per-pair correlation (unchanged from the original): a read x CpG 0/1
# indicator matrix per bucket, `cor(mat, use="pairwise.complete.obs")` for
# every CpG pair at once, matched with the number of reads actually
# covering BOTH members of that pair (phi/Pearson on a thin 2x2 table is
# unstable, hence MIN_READS_PER_PAIR). Each CpG PAIR contributes exactly
# one phi value to its distance bin (not weighted by supporting-read
# count) - final per-bin correlation is the mean of those pair-level phi
# values, matching how this metric is computed on the buffy coat side so
# the two curves stay comparable by construction.

# Given a fixed MAX_DISTANCE_BP << BUCKET_SIZE_BP, a CpG at position p
# belongs to its own bucket's gathering window trivially, and ALSO to the
# previous bucket's gathering window iff p is within MAX_DISTANCE_BP of
# that bucket's own boundary (p %% BUCKET_SIZE_BP < MAX_DISTANCE_BP) - at
# most 2 buckets per CpG, never more, since one full bucket-width of halo
# on one side is always enough when distances are capped well under a
# bucket width.
assign_gather_buckets <- function(dt, bucket_size_bp, max_distance_bp) {
  if (max_distance_bp >= bucket_size_bp) {
    stop("max_distance_bp must be < bucket_size_bp for the core+halo gathering logic to be correct")
  }
  dt <- copy(dt)
  dt[, bucket_idx := pos %/% bucket_size_bp]
  dt[, pos_in_bucket := pos %% bucket_size_bp]

  extra <- dt[pos_in_bucket < max_distance_bp]
  extra[, bucket_idx := bucket_idx - 1L]

  gathered <- rbind(dt, extra)
  gathered[, pos_in_bucket := NULL]
  gathered[, core_start := bucket_idx * bucket_size_bp]
  gathered[, core_end := (bucket_idx + 1L) * bucket_size_bp]
  gathered[, bucket_id := paste0(chrom, ":", core_start)]
  gathered
}

# Logged BEFORE the expensive per-bucket step so a crash recurrence leaves
# concrete evidence of which bucket/how large, instead of a blind OOM.
# Oversized buckets are DROPPED (not truncated/subsampled) - almost
# certainly a repeat-mapping/low-mappability pileup (satellite DNA,
# centromere, rDNA), not real coverage; a truncated sample of a pileup
# that size wouldn't produce a meaningful phi estimate anyway. Two
# SEPARATE caps: read count and CpG count, because the actual OOM
# mechanism (the read x CpG matrix, and separately the n_cpgs x n_cpgs
# cor()/distance matrices) scales differently in each dimension.
apply_bucket_size_caps <- function(gathered, max_reads_per_bucket, max_cpgs_per_bucket) {
  read_counts <- gathered[, .(n_reads = uniqueN(read_id)), by = bucket_id]
  log_msg("Reads-per-bucket: median=", median(read_counts$n_reads),
          " mean=", round(mean(read_counts$n_reads), 1),
          " p99=", quantile(read_counts$n_reads, 0.99), " max=", max(read_counts$n_reads))
  oversized_reads <- read_counts[n_reads > max_reads_per_bucket]$bucket_id
  if (length(oversized_reads) > 0) {
    log_msg("WARNING: ", length(oversized_reads), " bucket(s) exceed max_reads_per_bucket=",
            max_reads_per_bucket, " - dropping (likely repeat-mapping pileup):")
    print(read_counts[bucket_id %in% oversized_reads])
    gathered <- gathered[!bucket_id %in% oversized_reads]
  }

  cpg_counts <- gathered[, .(n_cpgs = uniqueN(pos)), by = bucket_id]
  log_msg("CpGs-per-bucket: median=", median(cpg_counts$n_cpgs),
          " mean=", round(mean(cpg_counts$n_cpgs), 1),
          " p99=", quantile(cpg_counts$n_cpgs, 0.99), " max=", max(cpg_counts$n_cpgs))
  oversized_cpgs <- cpg_counts[n_cpgs > max_cpgs_per_bucket]$bucket_id
  if (length(oversized_cpgs) > 0) {
    log_msg("WARNING: ", length(oversized_cpgs), " bucket(s) exceed max_cpgs_per_bucket=",
            max_cpgs_per_bucket, " - dropping (likely repeat/low-mappability artifact):")
    print(cpg_counts[bucket_id %in% oversized_cpgs])
    gathered <- gathered[!bucket_id %in% oversized_cpgs]
  }
  gathered
}

# One bucket's contribution: every CpG pair's phi + distance from this
# bucket's gathered read set, restricted to pairs whose lower position
# falls in [core_start, core_end).
compute_bucket_pair_correlation <- function(pos_vec, read_id_vec, state_vec,
                                             core_start_val, core_end_val,
                                             bin_width_bp, max_distance_bp, min_reads_per_pair) {
  empty <- data.table(distance_bin = integer(0), n_pairs = integer(0), sum_phi = numeric(0))

  mat <- data.table::dcast(data.table(read_id = read_id_vec, pos = pos_vec, state = state_vec),
                            read_id ~ pos, value.var = "state", fun.aggregate = function(x) x[1])
  mat <- as.matrix(mat[, -1, with = FALSE])
  if (ncol(mat) < 2L) return(empty)

  cor_mat <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))
  present <- (!is.na(mat)) * 1
  n_reads_covering_both <- crossprod(present)

  positions <- as.numeric(colnames(mat))
  dist_mat <- abs(outer(positions, positions, "-"))
  lower_pos_mat <- outer(positions, positions, pmin)

  ut <- upper.tri(cor_mat)
  phi <- cor_mat[ut]
  n_reads <- n_reads_covering_both[ut]
  distance <- dist_mat[ut]
  lower_pos <- lower_pos_mat[ut]

  keep <- !is.na(phi) & n_reads >= min_reads_per_pair & distance <= max_distance_bp &
    lower_pos >= core_start_val & lower_pos < core_end_val
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) return(empty)

  bin <- as.integer((distance[keep] %/% bin_width_bp) * bin_width_bp)
  data.table(distance_bin = bin, phi = phi[keep])[, .(n_pairs = .N, sum_phi = sum(phi)), by = distance_bin]
}

# Main entry point. `long_dt` must have columns: read_id, chrom, pos, state (0/1).
compute_pair_correlation_by_distance <- function(long_dt, max_distance_bp, bin_width_bp = 10L,
                                                  bucket_size_bp = 20000L, min_reads_per_pair = 5L,
                                                  max_reads_per_bucket = 20000L, max_cpgs_per_bucket = 3000L,
                                                  blacklist_dt = NULL) {
  dt <- copy(long_dt)

  if (!is.null(blacklist_dt)) {
    n_before <- nrow(dt)
    bl <- copy(blacklist_dt)
    setkey(bl, chrom, start, end)
    dt[, row_idx := .I]
    pts <- dt[, .(row_idx, chrom, start = pos, end = pos)]
    over <- foverlaps(pts, bl, by.x = c("chrom", "start", "end"), type = "any", nomatch = NULL)
    dt <- dt[!row_idx %in% unique(over$row_idx)]
    dt[, row_idx := NULL]
    log_msg("Rows after blacklist filtering: ", formatC(nrow(dt), format = "d", big.mark = ","),
            " / ", formatC(n_before, format = "d", big.mark = ","))
  }

  log_msg("Assigning core+halo gathering buckets (bucket_size=", bucket_size_bp,
          "bp, max_distance=", max_distance_bp, "bp)...")
  gathered <- assign_gather_buckets(dt, bucket_size_bp, max_distance_bp)
  gathered <- apply_bucket_size_caps(gathered, max_reads_per_bucket, max_cpgs_per_bucket)

  log_msg("Computing per-bucket CpG-pair correlation (", uniqueN(gathered$bucket_id), " buckets)...")
  per_bucket <- gathered[, compute_bucket_pair_correlation(
    pos, read_id, state, core_start[1], core_end[1], bin_width_bp, max_distance_bp, min_reads_per_pair
  ), by = bucket_id]

  result <- per_bucket[, .(n_pairs = sum(n_pairs), sum_phi = sum(sum_phi)), by = distance_bin]
  setorder(result, distance_bin)
  result[, mean_phi := sum_phi / n_pairs]
  result[, bin_mid := distance_bin + bin_width_bp / 2]
  result[]
}
