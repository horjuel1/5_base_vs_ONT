#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
Usage:
Rscript build_fivebase_site_table_from_modkit.r <modkit_cpg_bed> <out_parquet> [min_coverage]

Builds a 5-base per-CpG-site table from a `modkit pileup --cpg
--combine-strands` bedMethyl (produced by
build_fivebase_modkit_pileup.sh_launch.sh directly from the
cram_to_modbam.sh-produced modBAM), as an alternative to
build_fivebase_site_table.r's mod_pileup.txt-based route - both produce
the same output schema (chrom, cpg_pos, mod_count, coverage, beta), so
either can be pointed at 3_methylation_concordance/
compute_methylation_concordance.r without any changes there. Written to
compare against the rastair/mod_pileup.txt route: if this modkit-derived
table has substantially broader >=5x coverage, rastair call's
--min-baseq/--min-mapq thresholds (5CmC's downstream.sh) - not real
sequencing depth - were the bottleneck behind fivebase_sites.parquet's
thin site count.

Input is modkit's native (unreduced) bedMethyl - tab-separated, no
header, 18 columns (chrom, start, end, mod_code, score, strand,
thickStart, thickEnd, color, Nvalid_cov, percent_modified, Nmod,
Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall). --combine-strands
means modkit has already merged +/- strand calls and anchored each row at
the CpG's plus-strand C position - no manual minus-strand shift needed
here (unlike build_fivebase_site_table.r/build_ont_site_table.r, which
build cpg_pos from raw per-strand data). beta is recomputed as
Nmod/Nvalid_cov directly rather than trusting the reported
percent_modified column's scale/rounding, matching how beta is computed
everywhere else in this project.
")
}

bed_file     <- args[[1]]
out_file     <- args[[2]]
min_coverage <- if (length(args) >= 3) as.integer(args[[3]]) else 5L

if (!file.exists(bed_file)) stop("modkit_cpg_bed not found: ", bed_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", bed_file)
dt <- fread(bed_file, header = FALSE, sep = "\t",
            select = c(1, 2, 4, 10, 12),
            col.names = c("chrom", "cpg_pos", "mod_code", "coverage", "mod_count"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt <- dt[chrom %in% CHROM_ORDER & mod_code == "m"]
log_msg("5mC rows on autosomes: ", formatC(nrow(dt), format = "d", big.mark = ","))

site_dt <- dt[coverage >= min_coverage, .(chrom, cpg_pos, mod_count, coverage)]
site_dt[, beta := mod_count / coverage]
log_msg("CpG sites with coverage >= ", min_coverage, ": ", formatC(nrow(site_dt), format = "d", big.mark = ","))

write_parquet(site_dt, out_file)
log_msg("Saved: ", out_file)
