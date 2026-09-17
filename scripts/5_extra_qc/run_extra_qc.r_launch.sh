#!/bin/bash
# Depends on 2_persite_pileup/ outputs for coverage_breadth_overlap.r.
set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/5base_vs_ont"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$LOG_DIR"
module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_frags ont_pileup; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  PLOT_DIR="$OUT_BASE/$sample_id/plots"
  mkdir -p "$OUT_DIR" "$PLOT_DIR"

  sbatch \
    --job-name="extraqc_${sample_id}" \
    --nodes=1 --cpus-per-task=2 --mem=32G --time=4:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/extraqc_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/5_extra_qc/coverage_breadth_overlap.r' '$OUT_DIR/fivebase_sites.parquet' '$OUT_DIR/ont_sites.parquet' '$sample_id' '$OUT_DIR/coverage_breadth_overlap.tsv' && \
              Rscript '$PROJECT_DIR/scripts/5_extra_qc/snp_genotype_qc.r' '$fivebase_pileup' '$sample_id' '$OUT_DIR/snp_genotype_qc.tsv' && \
              Rscript '$PROJECT_DIR/scripts/5_extra_qc/methylation_by_fragment_length.r' '$fivebase_frag' '$ont_frags' '$sample_id' '$OUT_DIR/methylation_by_fragment_length.tsv' '$PLOT_DIR/methylation_by_fragment_length.png'"
done
