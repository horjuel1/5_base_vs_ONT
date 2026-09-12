#!/bin/bash
# Loops over config/samples.tsv and submits one sbatch job per sample.
set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/5base_vs_ont"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$OUT_BASE" "$LOG_DIR"
module load conda_R/4.4.x

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_fraglen ont_cpg_density ont_calls; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  mkdir -p "$OUT_DIR"

  sbatch \
    --job-name="ont_molbed_${sample_id}" \
    --nodes=1 --cpus-per-task=4 --mem=96G --time=1-00:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/ont_molbed_${sample_id}.o%j.txt" \
    --wrap="Rscript '$PROJECT_DIR/scripts/6_smmethid_compat/build_ont_moleculebed.r' '$ont_calls' '$OUT_DIR/${sample_id}_ont.molbed'"
done
