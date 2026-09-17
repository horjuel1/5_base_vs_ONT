#!/bin/bash
#SBATCH --job-name=rastair_recall
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --time=1-00:00:00
#SBATCH --partition=cancergen,shared
#SBATCH -o /dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison/logs/%x.o%j.txt

# Loops over config/samples.tsv, calling call_methylation_from_cram.sh per
# sample with a freshly-invoked rastair-2.1.1 (vs. whatever version 5CmC's
# own pipeline originally used - it pins v0.8.1). dedup.cram and
# mbias_masks.txt are derived from each sample's fivebase_frag directory,
# same convention as
# 6_smmethid_compat/build_fivebase_moleculebed.sh.

set -eo pipefail

PROJECT_DIR="/dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison"
FASTA="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/bwa/hg19_wm_idt_spike/hg19_wm_idt_spike.fa"
SAMPLES_TSV="$PROJECT_DIR/config/samples.tsv"
OUT_BASE="$PROJECT_DIR/results"

mkdir -p "$OUT_BASE"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id fivebase_frag fivebase_pileup ont_frags ont_pileup; do
  [[ -z "$sample_id" ]] && continue
  SAMPLE_DIR="$(dirname "$fivebase_frag")"
  CRAM="${SAMPLE_DIR}/dedup.cram"
  MASKS="${SAMPLE_DIR}/mbias_masks.txt"
  OUT_DIR="$OUT_BASE/$sample_id/tables/rastair_2.1.1_recall"

  bash "$(dirname "${BASH_SOURCE[0]}")/call_methylation_from_cram.sh" \
    "$CRAM" "$FASTA" "$OUT_DIR" "$MASKS" "${SLURM_CPUS_PER_TASK:-4}"
done
