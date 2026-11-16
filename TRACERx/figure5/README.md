# Figure 5 — SNV-phased allele-specific DMRs and intra-tumour methylation similarity

Split into **data generation (compute)** and **figure panels (plot)**. All plotting scripts are read-only on
`~/rdeocd/txMeth` and write to `outputs/`. Run everything from this directory.

| Panel | Content | Script |
|-------|---------|--------|
| **5A** | SNV-phased DMR classification pies (+ legend) | `scripts/plot_5A_SNV_phased_DMR_pies.R` |
| **5B–D** | Example SNV-phased DMRs (CDKN2A, MDGA1, CDH10) | `scripts/plot_5BCD_SNV_phased_DMR_examples.R` |
| **5E** | Bulk / purified tumour / normal methylation correlation heatmap (CRUK0062) | `scripts/plot_5E_region_correlation_heatmap.R` |
| **5F** | NJ sample trees (DMP m_b, DMP m_t, SNV CCF) | *(not in repo — made by Carla)* |
| **5G** | Δ(DMP − SNV) relative trunk length vs SBS4 | *(not in repo — made by Carla)* |

---

## 1. Data generation (phased methylation / SNV-phased DMRs)

The original PhD **notebook** is `~/data/CAMDAC_paper/figure5/scripts/Notes_phase_methylation_SNV_v2.R`
(not a CLI; paths such as `~/rdeocd/txMeth` are hardcoded; run interactively, section by section).
Only stages 1–3 below produce data; stage 4 is the original plotting code, replaced by the `plot_*` scripts here.

| Stage | Notebook lines | What it does | Output |
|-------|----------------|--------------|--------|
| 1. Phase methylation at SNVs | L1–437 | Per patient/region: phase CpG methylation to the mutant / wild-type allele at exome SNVs (`tx_SNV_calls.txt`, `GL_VAF < 0.1`), purified `m_t,1`/`m_t,2`, allele-specific DMP calls (`is_DM_t12`) | `txMeth/<patient>/SNVs/<R>/allele_specific_methylation_rates_<R>.RData`, `txMeth/SNV_DMPs.RData` |
| 2. SNV-phased DMR calling | L440–827 | Bin CpGs by annotation cluster (`hg19_all_regions_annotations.RData`), count (AS-)DMPs, consecutive DMPs, `DMR`, `DMR_1`, `DMR_2`, `CpG_phased_counts` per bin | `txMeth/<patient>/SNVs/<R>/SNV_phased_DMRs_calls_<R>.RData` (`df_phased_methylation`, `dt_anno_bins`) |
| 3. Cohort table | L828–881 | Collapse all samples' bins into one table (`m_t_1`, `m_t_2`, `m_*_DMPs`, `m_*_AS_DMPs`, `CN_1`, `CN_2`, …) | `txMeth/SNV_DMRs_v2.RData` (`dt`) |
| 4. Plots (superseded) | L884–1450 | Pies (L885–1146), enrichment tests, expression, example DMR plots (L1296–1450) | `txMeth_analysis/*.pdf` |

Inputs: `txMeth/tx_SNV_calls.txt`, `txMeth/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt`,
per-region CAMDAC results under `txMeth/<patient>/Methylation/<R>/`, `pipeline_files/loci_files/loci_hg19.RData`,
`pipeline_files/hg19_annotations/hg19_all_regions_annotations.RData`, `txMeth/tx421_drivers/*.csv`.

A related, cleaner SNV phasing script for the validation figure is `../figure2P-S_S6F-I/scripts/Phase_methylation_at_SNVs.R`
(this notebook's stage 1 is the DMP-oriented version).

Only the outputs of stage 2–3 are read by the plotting scripts below; they already exist for the TRACERx cohort
(`SNV_DMRs_v2.RData`, `SNV_phased_DMRs_calls_<R>.RData`), so stages 1–3 do **not** need to be rerun to remake the figure.

---

## 2. Plot panels

Shared helper: `scripts/fig5_common.R` (`classify_snv_dmrs()`, colour palette `snv_dmr_cols`).
Requires R packages `data.table`, `ggplot2`, `cowplot`, `optparse` (5E also `scales`).

### 5A–D: SNV-phased DMRs

Classification of a bin (`flag`): *trans* / *cis* / *LOH (loss of WT)* / *bi-allelic (both)* × hyper / hypo / none, from
`DMR_AS`, `DMR_1`, `DMR_2`, `CN_1`, `m_t_*_DMPs`, `m_t_*_AS_DMPs` versus `m_n_DMPs`. A DMR needs ≥ 3 consecutive
(AS-)DMPs; only bins with ≥ 3 phased CpGs (`CpG_phased_counts`) are used.

```bash
# 5A: pies + legend (counts: outputs/Fig5A_SNV_DMR_phasing_counts.tsv; 4,406 bins, 610 DMRs)
Rscript scripts/plot_5A_SNV_phased_DMR_pies.R \
  --snv_dmrs_rdata ~/rdeocd/txMeth/SNV_DMRs_v2.RData

# 5B–D: example DMRs (default: CDKN2A CRUK0082-R2, MDGA1 CRUK0071-R7, CDH10 CRUK0072-R4)
Rscript scripts/plot_5BCD_SNV_phased_DMR_examples.R \
  --tx_dir ~/rdeocd/txMeth \
  --snv_dmrs_rdata ~/rdeocd/txMeth/SNV_DMRs_v2.RData

# other loci: gene:patient:region[:class],...
Rscript scripts/plot_5BCD_SNV_phased_DMR_examples.R --examples "CDKN2A:CRUK0082:R2,MDGA1:CRUK0071:R7,CDH10:CRUK0072:R4:LOH_hypo"
```

CDH10 (5D) is a loss-of-WT bin whose cached `m_t_2_DMPs` is `NA`, so its class (`LOH_hypo`) is given explicitly.
The legend and cartoons in 5A had no code in the notebook; they are drawn in `plot_5A_SNV_phased_DMR_pies.R`.

Outputs: `outputs/Fig5A_SNV_phased_DMR_pies.{pdf,png}`, `outputs/Fig5A_SNV_DMR_phasing_counts.tsv`,
`outputs/Fig5BCD_SNV_phased_DMR_examples.{pdf,png}` and per-panel `outputs/Fig5{B,C,D}_SNV_phased_DMR_*.{pdf,png}`.

### 5E: correlation of bulk, purified tumour and normal methylation (CRUK0062)

Pairwise Pearson correlation between regions R1–R7 and the matched normal at all tumour-normal DMPs found in at least
one region (bulk DMPs for `m_b`, upper triangle; purified DMPs for `m_t`, lower triangle; normal in the top row and
right column). Samples are ordered by purity. Pass `--patient CRUK0062`; the script looks up the internal Patient
folder name from the tx100 clinical table.

```bash
Rscript scripts/plot_5E_region_correlation_heatmap.R \
  --tx_dir ~/rdeocd/txMeth \
  --patient CRUK0062

# check against the 2019 saved correlation table (max |difference| ~1e-15)
Rscript scripts/plot_5E_region_correlation_heatmap.R --recompute TRUE \
  --qc_workspace ~/rdeocd/txMeth/<Patient>/Clustering/clustering_workspace_<Patient>.RData
```

By default the script reads the per-region tables saved by the published analysis
(`txMeth/<Patient>/Clustering/clustering_workspace_<Patient>.RData`, `list_patient`), which reproduces the published
numbers. The current `Methylation/<R>/purified_tumour_DMP_calls.RData` files have been re-run since (different `m_t_corr`
and DMP sets; `m_t`-vs-normal correlations come out ~0.3 higher). Use `--clustering_workspace none` to recompute from
them instead. Method follows `R_scripts/Notes_clustering_DMPs_MVPs_pilot.R` (L19–165).
Correlations are cached in `outputs/Fig5E_region_correlations.tsv`; pass `--recompute TRUE` to rebuild.

Outputs: `outputs/Fig5E_region_correlation_heatmap.{pdf,png}`, `outputs/Fig5E_region_correlations.tsv`.

---

## Notes

- `~/rdeocd/txMeth/ASM_validation_Jul2021_v6.RData` and the other cached `txMeth` objects are never overwritten.
