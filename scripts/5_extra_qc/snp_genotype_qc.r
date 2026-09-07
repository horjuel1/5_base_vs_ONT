#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
source(file.path(script_dir, "..", "_shared_utils.r"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript snp_genotype_qc.r <mod_pileup_file> <sample_id> <out_file>

Rate of CpGs flagged as SNP-disrupted in the 5-base pileup (genotype !=
C/C on + strand, or != G/G on - strand) - the filter
2_persite_pileup/build_fivebase_site_table.r applies by default. A high
rate here would mean a large fraction of the genome is being silently
dropped before methylation concordance is even computed, which is worth
knowing rather than discovering implicitly from a smaller-than-expected
joined site count.
")
}

pileup_file <- args[[1]]
sample_id   <- args[[2]]
out_file    <- args[[3]]

dt <- fread(pileup_file, sep = "\t", select = c("#chr", "strand", "genotype"),
            col.names = c("chrom", "strand", "genotype"))
dt <- dt[chrom %in% CHROM_ORDER]

dt[, is_snp := !((strand == "+" & genotype == "C/C") | (strand == "-" & genotype == "G/G"))]

result <- data.table(
  sample_id = sample_id,
  n_positions = nrow(dt),
  n_flagged_snp = sum(dt$is_snp),
  pct_flagged_snp = 100 * mean(dt$is_snp)
)

fwrite_default(result, out_file)
print(result)
