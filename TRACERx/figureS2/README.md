# Supplementary Figure S2A–D

RRBS tumour LogR biases (MspI fragment length, GC content, replication timing)
and the combined natural-spline correction used by ASCAT.m.

| Panel | Content |
|-------|---------|
| A | Predicted LogR bias vs log10(MspI fragment length) |
| B | Predicted LogR bias vs fragment GC content |
| C | Predicted LogR bias vs replication timing |
| D | Combined correction: raw vs corrected LogR densities |

S2E–G (BAF rule schematics) are Illustrator / methods figures, not this script.

**Sources**
- Plotting: `R_scripts/archive/Notes_correcting_logR_Aug2018_Mar2019.R` (2019-04-02 Notes)
- Model: `LogR_correction()` in `R_scripts/run_ASCAT_RRBS_cluster_v3.R`
  (`lm(y ~ ns(length)+ns(GC)+ns(replic)`, `df=5`)

Uses per-sample `Copy_number/<sample>/<patient>.<sample>.SNPs.txt` (written after
LogR correction). Replication timing is re-annotated from
`pipeline_files/1000genomes_2012_v3_repliTiming/*.fst` (best-correlating track).

## Run (defaults = CRUK0062 / CRUK0062 R3)

```bash
Rscript scripts/plot_S2_LogR_bias_correction.R \
  --work_dir ~/data/CAMDAC_paper/figureS2 \
  --tx_dir ~/rdeocd/txMeth \
  --pipeline_dir ~/rdeocd/pipeline_files \
  --patient_id CRUK0062 \
  --sample_id R3 \
  --panels A,B,C,D
```

(`--patient_id` prefer CRUKID; Patient IDs also accepted; outputs are `CRUK0062_R3_S2*.pdf`.)

Requires the `fst` R package (for panel C/D replication annotation).

---

## Cohort: within-segment LogR noise vs correction

LogR correction does **not** change BAF. The relevant effect is **within-segment
LogR dispersion** on fixed ASCAT.m segments.

`scripts/summarise_LogR_noise_vs_correction.R`:

1. Prefer **resegmented** fragments:
   `Copy_number/<sample>/re_segmented/<patient>.<sample>.ascat.frag.resegmented.RData`
2. Single-region patients **CRUK0010** / **CRUK0047** (no reseg): use
   `Copy_number/<sample>/raw/<patient>.<sample>.ascat.frag.RData`
3. SNP-level raw / corrected LogR from `SNPs.txt` or `SNPs.RData`
4. Per autosomal segment (≥10 SNPs): `sd(LogR_raw)` vs `sd(LogR_corr)`;
   `Δ = noise_raw − noise_corr`
5. Length-weighted mean Δ per sample → mean ± SE per patient
6. Flag the **10 largest** segment Δ values cohort-wide

No ASPCF re-run; segment boundaries stay fixed (ideal for measuring LogR cleanup).

### Full cohort

```bash
Rscript scripts/summarise_LogR_noise_vs_correction.R \
  --work_dir ~/data/CAMDAC_paper/figureS2 \
  --tx_dir ~/rdeocd/txMeth
```

### Smoke test

```bash
Rscript scripts/summarise_LogR_noise_vs_correction.R \
  --work_dir ~/data/CAMDAC_paper/figureS2 \
  --tx_dir ~/rdeocd/txMeth \
  --patient_id CRUK0062,CRUK0047 \
  --sample_id R3,R2
```

### Outputs (`outputs/`)

| File | Content |
|------|---------|
| `LogR_noise_vs_correction_by_segment.tsv` | all segments (+ `top10` / `rank_delta`) |
| `LogR_noise_top10_segment_deltas.tsv` | 10 largest Δ segments |
| `LogR_noise_vs_correction_by_sample.tsv` | sample means |
| `LogR_noise_vs_correction_by_patient.tsv` | patient means ± SE |
| `LogR_noise_paired_by_sample.pdf` | raw vs corr (all samples; top-Δ labelled CRUKID_Rx) |
| `LogR_noise_paired_by_patient.pdf` | raw vs corr (patient means) |
| `LogR_noise_delta_by_patient_waterfall.pdf` | patient-averaged Δ (CRUKID on x) |

Patient / sample labels and TSV columns use **CRUKID** only (from
`txMeth/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt`;
override with `--clinical_info`). Internal TRACERx folder names use Patient IDs for
path discovery only.

### Note on earlier BAF script

`scripts/summarise_BAF_noise_vs_LogR_correction.R` re-ran ASPCF for BAF noise;
that is **not** the intended metric here (BAF should be unchanged under fixed
segments). Prefer the LogR script above.
