#!/bin/bash
#SBATCH --job-name=modkit_pileup_cpg
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=40G
#SBATCH --time=04-00:00:00
#SBATCH --partition=cancergen
#SBATCH -o /dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison/logs/%x.o%j.txt

# Regenerates a 5-base per-CpG pileup directly from the cram_to_modbam.sh
# -produced modBAM, using modkit's own --cpg pileup instead of
# either rastair's mod_pileup.txt (strict --min-baseq 30 --min-mapq 20 in
# 5CmC's downstream.sh - possibly why fivebase_sites.parquet's breadth
# looked thin) or bam2bed_5mC_5hmC_separate.sh's already-existing
# {sample}_5mC_all.bed (missing the strand column after its awk reduction,
# so it can't be correctly CpG-merged post-hoc without regenerating).
#
# --cpg WITHOUT --combine-strands: modkit's internal strand-combining
# panics on this modBAM (Rust panic in base_mods_adapter.rs, isolated to
# exactly the "combining strands at CpG motifs" step - basic MM/ML
# parsing, read sampling, and threshold estimation all worked fine before
# that point). Rather than chase a modkit bug, this leaves +/- strand
# rows separate; a parser doing the same manual minus-strand shift
# build_fivebase_site_table.r/build_ont_site_table.r already use would
# need to be paired with this (see
# 2_persite_pileup/build_ont_site_table_from_modkit.r for that exact
# pattern, currently wired to the ONT bed instead - this 5-base/modBAM
# path is not part of the current plan per the decision to keep 5-base
# on mod_pileup.txt, but the fix here is the same one that worked there).

set -eo pipefail

MODKIT="/dcs10/scharpf/data/horjuela/resources/modkit_v0.6.1/modkit"
MODBAM="/dcs11/scharpf/data/horjuela/sandbox/7_27_26_cram_to_modbam_LUCAS_full/modbams/CGPLLU431P_5B_lib1_bt2-bwa_umi_true_hg19_wm_idt_spike_v0.1.bam"
FASTA="/dcl01/scharpf/data/pipeline-hub/pipeline-resources/bwa/hg19_wm_idt_spike/hg19_wm_idt_spike.fa"
OUT_BED="/dcs11/scharpf/data/horjuela/9_2026_5-base_ONT_comparison/results/CGPLLU431P/tables/CGPLLU431P_5mC_per_strand.bed"

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
  --modified-bases 5mC \
  --threads "${SLURM_CPUS_PER_TASK:-4}"

echo "Wrote $(wc -l < "$OUT_BED") per-strand CpG rows to $OUT_BED"
