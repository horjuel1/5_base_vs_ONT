#!/bin/bash
# Loops over config/samples.tsv, submitting two sbatch jobs per sample
# (5-base pileup collapse is cheap; ONT calls.tsv aggregation is the
# heaviest single step in this pipeline - it's a per-call table, easily
# the largest raw file here, hence the larger mem/time request).
set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"
MIN_COVERAGE=5

mkdir -p "$OUT_BASE" "$LOG_DIR"

module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_fraglen ont_cpg_density ont_calls ont_pileup; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  mkdir -p "$OUT_DIR"

  sbatch \
    --job-name="fivebase_sites_${sample_id}" \
    --nodes=1 --cpus-per-task=2 --mem=32G --time=4:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/fivebase_sites_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/2_persite_pileup/build_fivebase_site_table.r' '$fivebase_pileup' '$OUT_DIR/fivebase_sites.parquet' $MIN_COVERAGE"

  sbatch \
    --job-name="ont_sites_${sample_id}" \
    --nodes=1 --cpus-per-task=4 --mem=96G --time=1-00:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/ont_sites_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/2_persite_pileup/build_ont_site_table.r' '$ont_pileup' '$OUT_DIR/ont_sites.parquet' $MIN_COVERAGE"
done
