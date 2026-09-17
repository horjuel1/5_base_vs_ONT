#!/bin/bash
# Loops over config/samples.tsv and submits one sbatch job per sample.
# Adjust PROJECT_DIR to wherever this repo actually lives on the cluster.
set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison"
SCRIPT="$PROJECT_DIR/scripts/1_fragment_lengths/compute_fragment_length_distribution.r"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$OUT_BASE" "$LOG_DIR"

module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_frags ont_pileup; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  mkdir -p "$OUT_DIR"

  sbatch \
    --job-name="fraglen_${sample_id}" \
    --nodes=1 \
    --cpus-per-task=2 \
    --mem=16G \
    --time=2:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/fraglen_${sample_id}.o%j.txt" \
    --wrap="Rscript '$SCRIPT' '$fivebase_frag' '$ont_frags' '$sample_id' '$OUT_DIR'"
done
