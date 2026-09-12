#!/bin/bash
# Loops over config/samples.tsv and submits one sbatch job per sample.
# Requires smMethID checked out on the cluster - set SMMETHID_DIR below.
# All current 5-base LUCAS samples share the same hg19_wm_idt_spike
# reference (per the sample directories' own naming); if/when a sample
# needs a different reference, move FASTA into a samples.tsv column
# instead of hardcoding it here.
set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/5base_vs_ont"
SMMETHID_DIR="/dcs11/scharpf/data/horjuela/smMethID"
FASTA="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/bwa/hg19_wm_idt_spike/hg19_wm_idt_spike.fa"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$OUT_BASE" "$LOG_DIR"
module load conda_R/4.4.x samtools

export SMMETHID_DIR

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_fraglen ont_cpg_density ont_calls; do
  [[ -z "$sample_id" ]] && continue
  OUT_DIR="$OUT_BASE/$sample_id/tables"
  mkdir -p "$OUT_DIR"

  sbatch \
    --job-name="fivebase_molbed_${sample_id}" \
    --nodes=1 --cpus-per-task=2 --mem=64G --time=1-00:00:00 \
    --partition=cancergen,shared \
    --output="$LOG_DIR/fivebase_molbed_${sample_id}.o%j.txt" \
    --export=ALL,SMMETHID_DIR="$SMMETHID_DIR" \
    --wrap="bash '$PROJECT_DIR/scripts/6_smmethid_compat/build_fivebase_moleculebed.sh' '$fivebase_frag' '$FASTA' '$OUT_DIR/${sample_id}_5base.molbed' '$OUT_DIR/${sample_id}_5base_frag.tsv' '$sample_id'"
done
