#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
Usage:
Rscript build_ont_moleculebed.r <ont_frags_file> <out_moleculebed_file> [min_call_prob]

Emits ONT's per-call table as a molecule-bed (chrom, start, end, read_id,
call [0/1]), header-less TSV - the same schema smMethID's own
bin/calls_to_moleculebed.r produces from a `modkit extract calls` TSV, and
the same schema its bin/mod_frags_to_moleculebed.r produces for 5-base
(see 6_smmethid_compat/build_fivebase_moleculebed.sh in this project).
One row per (read, CpG position) - not aggregated - so it can feed
smMethID's bin/molecule_bed_to_fragment_methylation.r directly, or be
row-bound with the 5-base molecule-bed for any tool that reads this
shared format.

Filtering/classification (secondary/supplementary-alignment exclusion,
autosome/within-alignment/fail filtering, m/h/- call_code classification,
minus-strand CpG anchoring) is identical to
2_persite_pileup/build_ont_site_table.r - both call the same
filter_and_classify_ont_calls() in _shared_utils.r, so a molecule-bed
built here and a site-pileup built there always agree on which calls
count.
")
}

calls_file    <- args[[1]]
out_file      <- args[[2]]
min_call_prob <- if (length(args) >= 3) as.numeric(args[[3]]) else NA_real_

if (!file.exists(calls_file)) stop("ont_frags_file not found: ", calls_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", calls_file)
dt <- fread(calls_file, sep = "\t",
            select = c("read_id", "chrom", "ref_position", "ref_mod_strand", "call_prob",
                       "call_code", "fail", "within_alignment", "flag"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt <- filter_and_classify_ont_calls(dt, min_call_prob = min_call_prob)

# A well-formed molecule-bed has one row per (read_id, cpg_pos) - the
# secondary/supplementary filter already removes the main source of
# duplicate/conflicting calls at the same position for one read_id, but
# collapse defensively (ANY_VALUE-style: keep the first) rather than
# silently emitting duplicate rows smMethID's own reducer doesn't expect.
molecule_bed <- unique(dt[, .(chrom, start = cpg_pos, end = cpg_pos + 1L, read_id, call = is_mod)],
                        by = c("chrom", "start", "read_id"))
setorder(molecule_bed, chrom, start)

fwrite(molecule_bed, out_file, sep = "\t", col.names = FALSE)
log_msg(sprintf("Wrote %s calls (%s unique molecules) to %s",
                 formatC(nrow(molecule_bed), format = "d", big.mark = ","),
                 formatC(uniqueN(molecule_bed$read_id), format = "d", big.mark = ","),
                 out_file))
