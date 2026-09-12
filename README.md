# 5-base vs ONT

Compares fragment length distribution, methylation concordance, longer-fragment
abundance, and CpG-to-CpG correlation by distance between 5-base (Illumina TAPS
short-read direct methylation) and ONT (Oxford Nanopore long-read) sequencing
of the same sample(s). Written in R, following the conventions of the sibling
`buffy_methylation_architecture` project (data.table/arrow/ggplot2, numbered
topic folders with `_launch.sh` SLURM wrappers, `commandArgs`-based scripts).

## Inputs (per sample)

| File | Platform | Description |
|---|---|---|
| `mod_frags.txt` | 5-base | per-fragment: insert size, per-CpG mod/unmod calls (offset-encoded) |
| `mod_pileup.txt` | 5-base | per-(position,strand) pileup: beta, mod/unmod counts, genotype |
| `*_fragment_lengths.tsv.gz` | ONT | per-read: fragment length, read length |
| `*_fragment_length_cpg_density.tsv.gz` | ONT | per-read: genomic/called/methylated CpG counts |
| `*_calls.tsv` | ONT | per-CpG-call-per-read: call probability, call code, alignment info |

Paths for known samples live in `config/samples.tsv`.

**Schema caveat**: the column schemas here were transcribed from pasted file
headers, not verified against the real files directly (this environment
can't reach the `/dcs05`/`/dcs11` cluster paths). Before a first real run,
sanity-check column count/order:
```
zcat -f <file> | head -1 | tr '\t' '\n' | cat -n
```
and adjust the `select=`/`col.names=` in the relevant script if anything
doesn't line up. `fail` in the ONT calls table is assumed to already flag
low-confidence calls to exclude — verify against whatever pipeline produced
`calls.tsv`.

## Pipeline stages

1. **`scripts/1_fragment_lengths/`** — fragment length distributions, mono-/
   di-nucleosome periodicity (reusing `buffy_methylation_architecture`'s
   already-verified peak-detection method), and longer-fragment abundance
   (including a DELFI-style short:long ratio, relevant since this is a LUCAS
   cfDNA sample).
2. **`scripts/2_persite_pileup/`** — collapses each platform's raw per-strand/
   per-call data into one row per CpG dinucleotide (chrom, cpg_pos, coverage,
   beta), so the platforms can be joined.
3. **`scripts/3_methylation_concordance/`** — joins the two per-site tables;
   reports Pearson/Spearman correlation, RMSD/MAD, Bland-Altman bias + limits
   of agreement, binary agreement/sensitivity/specificity/Cohen's kappa, and
   how these change with coverage depth.
4. **`scripts/4_cpg_correlation_distance/`** — within-read/fragment CpG
   co-methylation correlation vs. genomic distance, computed identically on
   both platforms via a shared core+halo bucketed engine (adapted from
   `buffy_methylation_architecture`'s `16_compute_cpg_pair_correlation_
   unconstrained_per_sample.r`). This is where ONT's long-read advantage
   should show up most clearly — short fragments/reads can't be measured
   past their own length.
5. **`scripts/5_extra_qc/`** — coverage breadth/overlap between platforms,
   5-base SNP/genotype QC (rate of CpGs dropped as SNP-disrupted), and
   methylation/CpG-density as a function of fragment length (relevant for
   cfDNA nucleosome-protection signal).
6. **`scripts/6_smmethid_compat/`** — interoperability with the sibling
   `smMethID` project's shared "molecule-bed" format (`chrom, start, end,
   read_id, call`), which both an ONT and a 5-base converter there already
   emit. `build_ont_moleculebed.r` emits it from our `calls.tsv` using the
   same filtering as stage 2. `build_fivebase_moleculebed.sh` wraps
   smMethID's own `mod_frags_to_moleculebed.r` (CIGAR-aware CpG position
   decoding, unlike stages 2/4's naive `start + offset`) and its
   `molecule_bed_to_fragment_methylation.r` reducer, deriving the sample's
   `dedup.cram`/`mbias_masks.txt` paths automatically from `mod_frags.txt`'s
   own directory (5CmC's per-sample layout). `mod_frags_to_moleculebed.r`
   now refuses to write output if fewer than `--min_validation_rate`
   (default 90%) of raw offsets resolve to a confirmed CG dinucleotide in
   the reference — watch that line in the log on a real sample; a low rate
   means the offset convention doesn't match your rastair version.

   Requires `SMMETHID_DIR` set to wherever smMethID lives on the cluster,
   `samtools`, and the R `argparser` package.

Run order: 1 and 2 have no dependencies on each other; 3 and 5's
`coverage_breadth_overlap.r` depend on 2; 4 is independent of 2/3; 6 is
independent of everything else (reads the same raw inputs directly).

## Running

Each `.r` script is self-documenting — run with no arguments to print usage.
Locally (no cluster access), e.g.:
```
Rscript scripts/1_fragment_lengths/compute_fragment_length_distribution.r \
  <mod_frag.txt> <fragment_lengths.tsv.gz> <sample_id> <out_dir>
```

On the cluster, each stage has a `*_launch.sh` that loops over
`config/samples.tsv` and submits one `sbatch` job per sample. **Set
`PROJECT_DIR` at the top of each `_launch.sh` to wherever this repo actually
lives on the cluster** before submitting — it's currently a placeholder
guess (`/dcs11/scharpf/data/horjuela/5base_vs_ont`).

Requires R with `data.table`, `arrow`, `ggplot2`, `patchwork` (all already
available via `module load conda_R/4.4.x` on the cluster, per the sibling
project). `PlasmaTools` (for hg19 blacklist filtering in stage 4's ONT
correlation) is optional — used if installed, skipped with a warning
otherwise. Stage 6's 5-base wrapper additionally needs `samtools`, the R
`argparser` package, and a `smMethID` checkout (path via `SMMETHID_DIR`).

## Known simplifications / things to verify against real data

- 5-base's `mod_cpgs`/`unmod_cpgs` offsets are assumed 0-based from `start`
  in stages 2 and 4 (both mates of a `read_id` pooled into one fragment-level
  CpG list); this is a naive approximation that is silently wrong across an
  indel or soft-clip (rastair's offsets are query-cycle positions, not plain
  reference offsets). For anything where exact 5-base CpG position matters,
  use `scripts/6_smmethid_compat/build_fivebase_moleculebed.sh` instead,
  which does full CIGAR-aware decoding via smMethID's converter — validated
  end-to-end here against a synthetic ref+CRAM+mod_frags.txt with known
  planted CpG positions (see that script's header comment), though not yet
  against a real cluster sample.
- ONT `call_code`'s first character is assumed `m`=5mC, `h`=5hmC,
  `-`=canonical; m+h are counted together as "modified" to match TAPS
  chemistry (which converts both 5mC and 5hmC). Calls from secondary/
  supplementary alignment records (SAM flag 0x100/0x800) are dropped before
  any per-read or per-site aggregation, matching smMethID's own
  `calls_to_moleculebed.r` — a chimeric/concatemer ONT read can otherwise
  pool CpG calls from two physically unrelated genomic loci under one
  `read_id`.
- Both platforms are assumed hg19-aligned (per the 5-base folder name) — no
  liftover is implemented; if ONT turns out to be hg38, site-level joins in
  stages 3 and 5 will silently produce near-zero overlap rather than erroring.
- Stage 4's core+halo bucketing requires `max_distance_bp < BUCKET_SIZE_BP`
  (20,000bp) — raise `BUCKET_SIZE_BP` in the two compute scripts if you need
  a wider distance range than that.
