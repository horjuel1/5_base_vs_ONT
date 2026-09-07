#!/bin/bash
# Depends on 2_persite_pileup/build_persite_pileups.r_launch.sh having
# already produced fivebase_sites.parquet and ont_sites.parquet per sample.
set -euo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/5base_vs_ont"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$LOG_DIR"
module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_fraglen ont_cpg_density ont_calls; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"

  sbatch \
    --job-name="concordance_${sample_id}" \
    --nodes=1 --cpus-per-task=2 --mem=32G --time=4:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/concordance_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/3_methylation_concordance/compute_methylation_concordance.r' '$OUT_DIR/fivebase_sites.parquet' '$OUT_DIR/ont_sites.parquet' '$sample_id' '$OUT_DIR'"
done
