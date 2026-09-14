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
Rscript build_fivebase_site_table_from_modkit.r <modkit_per_strand_cpg_bed> <out_parquet> [min_coverage]

Builds a 5-base per-CpG-site table from a `modkit pileup --cpg` bedMethyl
WITHOUT --combine-strands (produced by
build_fivebase_modkit_pileup.sh_launch.sh directly from the
cram_to_modbam.sh-produced modBAM) - as an alternative to
build_fivebase_site_table.r's mod_pileup.txt-based route - both produce
the same output schema (chrom, cpg_pos, mod_count, coverage, beta), so
either can be pointed at 3_methylation_concordance/
compute_methylation_concordance.r without any changes there. Written to
compare against the rastair/mod_pileup.txt route: if this modkit-derived
table has substantially broader >=5x coverage, rastair call's
--min-baseq/--min-mapq thresholds (5CmC's downstream.sh) - not real
sequencing depth - were the bottleneck behind fivebase_sites.parquet's
thin site count.

--combine-strands was dropped from the modkit invocation because it
panics (Rust 'entered unreachable code' in base_mods_adapter.rs) on this
particular modBAM's MM/ML tags - so this script does the same manual
minus-strand shift build_fivebase_site_table.r and
build_ont_site_table.r already do, instead of depending on modkit's
internal (here, crashing) strand-combining.

Input is modkit's native (unreduced) bedMethyl - tab-separated, no
header, 18 columns (chrom, start, end, mod_code, score, strand,
thickStart, thickEnd, color, Nvalid_cov, percent_modified, Nmod,
Ncanonical, Nother_mod, Ndelete, Nfail, Ndiff, Nnocall). Without
--combine-strands, a CpG's '+' strand row (the C) sits at `start`; its
'-' strand row (the G, i.e. the C measured on the opposite strand) sits
at `start+1` - same convention as mod_pileup.txt - so the '-' strand row
is shifted back 1bp to land on the same cpg_pos as the '+' strand row
before summing coverage/mod_count across both. beta is computed as
Nmod/Nvalid_cov directly (summed across both strands) rather than
trusting the reported percent_modified column's scale/rounding, matching
how beta is computed everywhere else in this project.
")
}

bed_file     <- args[[1]]
out_file     <- args[[2]]
min_coverage <- if (length(args) >= 3) as.integer(args[[3]]) else 5L

if (!file.exists(bed_file)) stop("modkit_per_strand_cpg_bed not found: ", bed_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", bed_file)
dt <- fread(bed_file, header = FALSE, sep = "\t",
            select = c(1, 2, 4, 6, 10, 12),
            col.names = c("chrom", "start", "mod_code", "strand", "coverage", "mod_count"))
log_msg("Rows read: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt <- dt[chrom %in% CHROM_ORDER & mod_code == "m"]
log_msg("5mC rows on autosomes: ", formatC(nrow(dt), format = "d", big.mark = ","))

dt[, cpg_pos := ifelse(strand == "-", start - 1L, start)]

site_dt <- dt[, .(mod_count = sum(mod_count), coverage = sum(coverage)), by = .(chrom, cpg_pos)]
site_dt <- site_dt[coverage >= min_coverage]
site_dt[, beta := mod_count / coverage]
log_msg("CpG sites with coverage >= ", min_coverage, ": ", formatC(nrow(site_dt), format = "d", big.mark = ","))

write_parquet(site_dt, out_file)
log_msg("Saved: ", out_file)
