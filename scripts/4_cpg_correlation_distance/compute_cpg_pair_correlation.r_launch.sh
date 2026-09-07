#!/bin/bash
# The ONT correlation step is the heaviest job in this whole pipeline (a
# per-bucket read x CpG matrix + cor() over the full calls.tsv) - matches
# the resource profile buffy_methylation_architecture's own
# 16_compute_cpg_pair_correlation_unconstrained_per_sample.r needed on
# real cluster data after its OOM history, hence the generous mem/time.
set -euo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/5base_vs_ont"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"
MAX_DISTANCE_BP=2000

mkdir -p "$LOG_DIR"
module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_fraglen ont_cpg_density ont_calls; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  PLOT_DIR="$OUT_BASE/$sample_id/plots"
  mkdir -p "$OUT_DIR" "$PLOT_DIR"

  sbatch \
    --job-name="cpgcorr_5base_${sample_id}" \
    --nodes=1 --cpus-per-task=2 --mem=32G --time=8:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/cpgcorr_5base_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/4_cpg_correlation_distance/1_compute_fivebase_cpg_pair_correlation.r' '$fivebase_frag' '$sample_id' '$OUT_DIR/cpg_correlation_5base.tsv' $MAX_DISTANCE_BP"

  sbatch \
    --job-name="cpgcorr_ont_${sample_id}" \
    --nodes=1 --cpus-per-task=4 --mem=96G --time=1-00:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/cpgcorr_ont_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/4_cpg_correlation_distance/2_compute_ont_cpg_pair_correlation.r' '$ont_calls' '$sample_id' '$OUT_DIR/cpg_correlation_ont.tsv' $MAX_DISTANCE_BP"
done

# Run 3_plot_correlation_by_distance_comparison.r manually per sample once
# both jobs above have finished (it's a quick plot, not worth its own sbatch).
