#!/bin/bash
#SBATCH --job-name=modkit_pileup_cpg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04-00:00:00
#SBATCH --partition=cancergen
#SBATCH -o /dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison/logs/%x.o%j.txt

# Regenerates a 5-base per-CpG pileup directly from the cram_to_modbam.sh
# -produced modBAM, using modkit's own --cpg --combine-strands instead of
# either rastair's mod_pileup.txt (strict --min-baseq 30 --min-mapq 20 in
# 5CmC's downstream.sh - possibly why fivebase_sites.parquet's breadth
# looked thin) or bam2bed_5mC_5hmC_separate.sh's already-existing
# {sample}_5mC_all.bed (missing the strand column after its awk reduction,
# so it can't be correctly CpG-merged post-hoc without regenerating).
#
# --cpg --combine-strands: modkit handles +/- strand merging internally
# and anchors each combined CpG at the plus-strand C's position - the
# same coordinate convention this project's build_ont_site_table.r and
# build_fivebase_site_table.r already use, so the output here needs no
# manual minus-strand shift before joining against ont_sites.parquet.

set -eo pipefail

MODKIT="/dcs10/scharpf/data/horjuela/resources/modkit_v0.6.1/modkit"
MODBAM="/dcs11/scharpf/data/horjuela/sandbox/7_27_26_cram_to_modbam_LUCAS_full/modbams/CGPLLU431P_5B_lib1_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1.bam"
FASTA="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/bwa/hg19_wm_idt_spike/hg19_wm_idt_spike.fa"
OUT_BED="/dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison/results/CGPLLU431P/tables/CGPLLU431P_5mC_cpg_combined.bed"

mkdir -p "$(dirname "$OUT_BED")"

for f in "$MODKIT" "$MODBAM" "$FASTA"; do
  [[ -e "$f" ]] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
[[ -f "${MODBAM}.bai" || -f "${MODBAM%.bam}.bai" ]] || echo "WARNING: no BAM index found for $MODBAM" >&2

"$MODKIT" pileup \
  "$MODBAM" \
  "$OUT_BED" \
  --ref "$FASTA" \
  --cpg \
  --combine-strands \
  --modified-bases 5mC \
  --threads "${SLURM_CPUS_PER_TASK:-4}"

echo "Wrote $(wc -l < "$OUT_BED") CpG sites to $OUT_BED"
