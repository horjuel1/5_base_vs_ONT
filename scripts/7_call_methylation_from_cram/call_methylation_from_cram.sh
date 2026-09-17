#!/bin/bash
# Calls methylation directly from a CRAM/BAM using rastair, producing both
# a CpG-site pileup (mod_pileup.txt-equivalent) and per-read methylation
# calls (mod_frags.txt-equivalent) - the same two outputs 5CmC's own
# pipeline_scripts/downstream.sh produces, but standalone/reusable here
# with a specific, controllable rastair binary rather than depending on a
# full 5CmC pipeline run having already happened for a sample.
#
# Written to test whether rastair-2.1.1 (available at
# /dcs04/scharpf/data/pipeline-hub/pipeline-lib/rastair-2.1.1, per a
# colleague's own methylation_5base.sh) gives materially different CpG
# coverage/breadth than whatever rastair version originally produced this
# sample's mod_pileup.txt/mod_frags.txt (5CmC pins v0.8.1) - part of the
# ongoing "why are there so few shared 5-base/ONT sites" investigation.
# Same rationale as the modkit-on-modBAM detour
# (2_persite_pileup/build_fivebase_modkit_pileup.sh_launch.sh), but
# testing the rastair version itself rather than a different tool/pipeline
# entirely.
#
# Output naming matches 5CmC's own convention exactly
# ({out_dir}/mod_pileup.txt, {out_dir}/mod_frags.txt) so
# build_fivebase_site_table.r, 1_compute_fivebase_cpg_pair_correlation.r,
# and build_fivebase_moleculebed.sh can all read this output directly, no
# changes needed.
#
# Usage:
#   call_methylation_from_cram.sh <cram_or_bam> <fasta_file> <out_dir> [mbias_masks_file] [threads]
#
# mbias_masks_file (optional): 5CmC-format CSV with "OT,L,R,L,R" / "OB,L,R,L,R"
# rows (same file build_fivebase_moleculebed.sh already reads) - only used
# for `rastair call`'s --nOT/--nOB, matching 5CmC's downstream.sh convention
# of masking the site-level pileup but NOT the per-read output (mod_frags.txt
# is documented there as "no masking" - masks are applied post-hoc during
# offset decoding instead, e.g. by mod_frags_to_moleculebed.r). Omit or pass
# a nonexistent path to skip masking (--nOT/--nOB 0,0,0,0).

set -eo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cram_or_bam> <fasta_file> <out_dir> [mbias_masks_file] [threads]" >&2
  exit 1
fi

BAM="$1"
FASTA="$2"
OUT_DIR="$3"
MASKS="${4:-}"
THREADS="${5:-4}"

RASTAIR="${RASTAIR:-/dcs04/scharpf/data/pipeline-hub/pipeline-lib/rastair-2.1.1}"

for f in "$BAM" "$FASTA" "$RASTAIR"; do
  [[ -e "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"

OT_MASK="0,0,0,0"
OB_MASK="0,0,0,0"
if [[ -n "$MASKS" && -f "$MASKS" ]]; then
  OT_MASK=$(awk -F',' '$1=="OT"{print $2","$3","$4","$5}' "$MASKS")
  OB_MASK=$(awk -F',' '$1=="OB"{print $2","$3","$4","$5}' "$MASKS")
  [[ -z "$OT_MASK" ]] && OT_MASK="0,0,0,0"
  [[ -z "$OB_MASK" ]] && OB_MASK="0,0,0,0"
  echo "Using masks from ${MASKS}: OT=${OT_MASK} OB=${OB_MASK}" >&2
else
  echo "No mbias_masks_file given/found - calling with no masking (OT=OB=0,0,0,0)." >&2
fi

echo "Step 1: rastair call -> ${OUT_DIR}/mod_pileup.txt" >&2
"$RASTAIR" call \
  --nOT "$OT_MASK" \
  --nOB "$OB_MASK" \
  --min-baseq 30 \
  --min-mapq 20 \
  --fasta-file "$FASTA" \
  -@ "$THREADS" \
  "$BAM" \
  > "${OUT_DIR}/mod_pileup.txt"

echo "Step 2: rastair per-read -> ${OUT_DIR}/mod_frags.txt" >&2
"$RASTAIR" per-read \
  --min-mapq 20 \
  --fasta-file "$FASTA" \
  -@ "$THREADS" \
  "$BAM" \
  > "${OUT_DIR}/mod_frags.txt"

echo "Done." >&2
echo "  $(wc -l < "${OUT_DIR}/mod_pileup.txt") pileup rows -> ${OUT_DIR}/mod_pileup.txt" >&2
echo "  $(wc -l < "${OUT_DIR}/mod_frags.txt") fragment rows -> ${OUT_DIR}/mod_frags.txt" >&2
