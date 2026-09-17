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
Rscript build_ont_site_table_from_modkit.r <ont_modkit_bed> <out_parquet> [min_coverage]

Builds the ONT per-CpG-site table from a modkit-pileup-derived bed
(config/samples.tsv's ont_pileup column), as an alternative to
build_ont_site_table.r's calls.tsv-based route - both produce the same
output schema (chrom, cpg_pos, mod_count, coverage, beta), so either can
be pointed at 3_methylation_concordance/compute_methylation_concordance.r
without any changes there.

Input is a REDUCED bedMethyl, tab-separated, no header, 8 columns in
this exact order: chrom, start, end, mod_code, strand, coverage
(Nvalid_cov), percent_modified, mod_count (Nmod). This is NOT modkit's
native 18-column layout and NOT the original lab convention's 7-column
reduction either (bam2bed_5mC_5hmC_separate.sh's own
'$1,$2,$3,$4,$12,$10,$11' omits strand entirely, which silently breaks
CpG-dinucleotide merging - the exact bug this schema exists to avoid).
Produce it with:
  awk 'BEGIN {OFS=\"\t\"} $4 == \"m\" {print $1,$2,$3,$4,$6,$10,$11,$12}' <raw_modkit_pileup_bed>
(raw_modkit_pileup_bed = modkit pileup's own unreduced output, i.e. skip
--combine-strands and the awk step that drops column 6/strand).

A CpG's '+' strand row (the C) sits at `start`; its '-' strand row (the
G, i.e. the C measured on the opposite strand) sits at `start+1` - same
convention as mod_pileup.txt and build_ont_site_table.r - so the '-'
strand row is shifted back 1bp to land on the same cpg_pos as the '+'
strand row before summing coverage/mod_count across both. beta is
computed as mod_count/coverage directly (summed across both strands)
rather than trusting the reported percent_modified column's scale/
rounding, matching how beta is computed everywhere else in this project.
")
}

bed_file     <- args[[1]]
out_file     <- args[[2]]
min_coverage <- if (length(args) >= 3) as.integer(args[[3]]) else 5L

if (!file.exists(bed_file)) stop("ont_modkit_bed not found: ", bed_file)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading: ", bed_file)
dt <- fread(bed_file, header = FALSE, sep = "\t",
            select = c(1, 2, 4, 5, 6, 8),
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
