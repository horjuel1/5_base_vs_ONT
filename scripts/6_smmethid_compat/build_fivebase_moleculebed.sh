#!/bin/bash
# Converts one 5-base sample's mod_frags.txt into a molecule-bed using
# smMethID's own CIGAR-aware converter (bin/mod_frags_to_moleculebed.r),
# instead of this project's own naive start+offset position decoding
# (2_persite_pileup/build_fivebase_site_table.r,
# 4_cpg_correlation_distance/1_compute_fivebase_cpg_pair_correlation.r).
# Those two scripts are unaffected/still run - this is an additional,
# more accurate path for anything that needs true per-CpG genomic
# position on the 5-base side (e.g. joining with smMethID's own
# molecule_bed_to_fragment_methylation.r, or eventually replacing this
# project's own offset math).
#
# VALIDATION: end-to-end tested against a hand-built ref.fa + CRAM +
# mod_frags.txt (planted CpGs at known positions, no indel/soft-clip) -
# mod_frags_to_moleculebed.r correctly recovered the exact planted 0-based
# genomic coordinates through this wrapper's full mask-derivation ->
# samtools-CRAM-read -> molecule-bed -> fragment-methylation chain. This
# is in addition to (not a replacement for) mod_frags_to_moleculebed.r's
# own toy CIGAR test suite (indels/soft-clips/insertions) - still NOT yet
# run against a real cluster sample's mod_frags.txt/BAM.
#
# On a REAL sample, watch the "CpG offset validation" line in this job's
# log: mod_frags_to_moleculebed.r now refuses to write output below
# --min_validation_rate (default 90%) of raw mod_cpgs/unmod_cpgs offsets
# resolving to a confirmed CG dinucleotide in the reference - a low rate
# means this sample's rastair version/column layout doesn't match what
# the script assumes (see smMethID/docs/method-improvement/
# 17-mod-frags-direct-input.md), and it stops rather than silently
# writing wrong positions. Verified this check fires correctly on a
# deliberately shifted (wrong) offset test case (0% validation, refused
# to write) as well as the correct one (100%, wrote normally).
#
# BAM and m-bias mask paths are derived from the mod_frags.txt file's own
# directory, matching 5CmC's per-sample output convention (same
# derivation smMethID/launch_mod_frags_to_moleculebed_array.sh uses):
#   <sample_dir>/dedup.cram        (source BAM/CRAM)
#   <sample_dir>/mbias_masks.txt   (OT/OB m-bias masks, CSV: "OT,L,R,L,R" / "OB,L,R,L,R")
#
# Usage:
#   build_fivebase_moleculebed.sh <mod_frags_file> <fasta_file> <out_moleculebed_file> <out_fragtsv_file> [sample_id]
#
# <fasta_file> must be the SAME reference the sample was aligned to - for
# this project's CGPLLU431P, that's the hg19_wm_idt_spike reference named
# in the sample's own directory
# (.../CGPLLU431P_5B_lib1_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1/), e.g.
# /dcl01/scharpf/data/pipeline-hub/pipeline-resources/bwa/hg19_wm_idt_spike/hg19_wm_idt_spike.fa
# per smMethID's own launch_mod_frags_to_moleculebed_array.sh usage example
# - confirm this path still resolves on the cluster before relying on it.
set -eo pipefail

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <mod_frags_file> <fasta_file> <out_moleculebed_file> <out_fragtsv_file> [sample_id] [min_validation_rate]" >&2
  exit 1
fi

FRAG_FILE="$1"
FASTA="$2"
MOLBED="$3"
FRAGTSV="$4"
SAMPLE_ID="${5:-}"
MIN_VALIDATION_RATE="${6:-0.9}"

SMMETHID_DIR="${SMMETHID_DIR:-/dcs11/scharpf/data/horjuela/smMethID}"
SAMPLE_DIR="$(dirname "$FRAG_FILE")"
BAM_FILE="${SAMPLE_DIR}/dedup.cram"
MASKS="${SAMPLE_DIR}/mbias_masks.txt"

for f in "$FRAG_FILE" "$FASTA" "$BAM_FILE"; do
  [[ -f "$f" ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done
[[ -f "${SMMETHID_DIR}/bin/mod_frags_to_moleculebed.r" ]] || {
  echo "ERROR: smMethID not found at SMMETHID_DIR=${SMMETHID_DIR} - set the env var to wherever smMethID lives on the cluster." >&2
  exit 1
}

# Same awk parse as smMethID/bin/cram_to_modbam.sh and
# launch_mod_frags_to_moleculebed_array.sh: "OT,L,R,L,R" / "OB,L,R,L,R"
# rows in mbias_masks.txt -> "L,R,L,R" for --ot_mask/--ob_mask. Falls back
# to no masking (with a warning) if the file isn't present, rather than
# failing outright - masking is an accuracy refinement, not a hard
# requirement for the conversion to run.
OT_MASK="0,0,0,0"
OB_MASK="0,0,0,0"
if [[ -f "$MASKS" ]]; then
  OT_MASK=$(awk -F',' '$1=="OT"{print $2","$3","$4","$5}' "$MASKS")
  OB_MASK=$(awk -F',' '$1=="OB"{print $2","$3","$4","$5}' "$MASKS")
  [[ -z "$OT_MASK" ]] && OT_MASK="0,0,0,0"
  [[ -z "$OB_MASK" ]] && OB_MASK="0,0,0,0"
  echo "Using masks from ${MASKS}: OT=${OT_MASK} OB=${OB_MASK}" >&2
else
  echo "WARNING: ${MASKS} not found - proceeding with no m-bias masking (OT=OB=0,0,0,0)." >&2
fi

mkdir -p "$(dirname "$MOLBED")" "$(dirname "$FRAGTSV")"

echo "Step 1: mod_frags.txt -> molecule-bed" >&2
Rscript "${SMMETHID_DIR}/bin/mod_frags_to_moleculebed.r" \
  --frag_file "$FRAG_FILE" \
  --fasta_file "$FASTA" \
  --bam_file "$BAM_FILE" \
  --output_file "$MOLBED" \
  --ot_mask "$OT_MASK" --ob_mask "$OB_MASK" \
  --min_validation_rate "$MIN_VALIDATION_RATE"

echo "Step 2: molecule-bed -> per-fragment methylation table" >&2
SAMPLE_ARG=(); [[ -n "$SAMPLE_ID" ]] && SAMPLE_ARG=(--sample_id "$SAMPLE_ID")
Rscript "${SMMETHID_DIR}/bin/molecule_bed_to_fragment_methylation.r" \
  --molecule_bed "$MOLBED" \
  --length_file "$FRAG_FILE" \
  --platform 5base "${SAMPLE_ARG[@]}" \
  --output_file "$FRAGTSV"

echo "Wrote ${MOLBED} and ${FRAGTSV}" >&2
