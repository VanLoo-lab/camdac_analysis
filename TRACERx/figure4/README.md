# Figure 4 — allele-specific methylation (ASM)

Pipeline split into **phasing / classification (compute)** vs **figure panels (plot)**.

## Compute (no figure output)

| Step | Script | Output |
|------|--------|--------|
| 1. Phase CpGs to `genotype_1` / `genotype_2` | `../figure2P-S_S6F-I/scripts/Phase_methylation_at_SNVs.R` | `allele_specific_methylation_rates_{sample}.RData` |
| 2. Power + clonal/subclonal DMP classes | `scripts/compute_DMP_classification.R` | `DMP_classification_{sample}.RData` |
| 3. Cohort DMP count tables (Fig 4I/J) | `scripts/compute_ASM_cohort_summary_figure4J.R` | `outputs/ASM_cohort_summary_figure4J_refresh.RData` |

```bash
# Per cohort (tumor regions R*; requires txMeth mount)
Rscript scripts/compute_DMP_classification.R \
  --tx_dir ~/rdeocd/txMeth \
  --work_dir ~/data/CAMDAC_paper/figure4 \
  --n_cores 8
```

Shared helpers: `scripts/asm_utils.R` (HDI, beta overlap, `h()` DMP test).

## Plot panels

| Panel | Content | Script | Original PhD source |
|-------|---------|--------|---------------------|
| **4A** | Imprinting pie charts by patient × locus | `scripts/plot_4A_imprinting.R` | `Notes_imprinted_genes_v2.R` (~L729–1025); data: `imprinting_summary_tumours_Oct2021_v3.RData` |
| **4B** | Phylogenetic trees + region status | *(manual figure)* | Hand-annotated; no R reproduction in repo |
| **4C–D** | Example region: CpG count histograms + `m_t,1` vs `m_t,2` | `scripts/plot_4CDEF_example_tumour.R` | `Notes_ASM_plots_14Jul2021.R` (~L191–291, per-sample `*_ASM_validation_v2.pdf`) |
| **4E–F** | Multi-region `m_t` and `m_t,1−m_t,2` comparisons | `scripts/plot_4CDEF_example_tumour.R` | `Notes_ASM_plots_14Jul2021.R` (~L755–924) |
| **4G** | SNV phylogeny + cancer cell fraction | *(manual figure)* | Hand-annotated; no R reproduction in repo |
| **4H** | UMAP coloured by epimutation class | `scripts/plot_4H_umap.R` | `Notes_ASM_plots_14Jul2021.R` (~L322–599); loads `outputs/ASM_UMAP_recovered_from_pdf.RData` |
| **4I** | Cohort DMP assignment bars | `scripts/plot_4IJ_cohort_DMP_summary.R --panels I` | `Notes_ASM_plots_14Jul2021.R` / `Notes_ASM_plots_16Jun2021.R` (~L987–1042) |
| **4J** | Clinical feature odds ratios | `scripts/plot_4IJ_cohort_DMP_summary.R --panels J` | `Notes_ASM_plots_14Jul2021.R` (~L1166–1224) |

### Regenerate extracted panels

```bash
# 4A (imprinting pie overview — not the S8C bar summary)
Rscript scripts/plot_4A_imprinting.R --tx_dir ~/rdeocd/txMeth

# 4H: load stored UMAP coordinates / classes (default CRUK0084)
Rscript scripts/plot_4H_umap.R --tx_dir ~/rdeocd/txMeth --crukid CRUK0084

# 4I–J: refresh cohort tables (does not overwrite ASM_validation_Jul2021_v6.RData)
Rscript scripts/compute_ASM_cohort_summary_figure4J.R --tx_dir ~/rdeocd/txMeth

Rscript scripts/plot_4IJ_cohort_DMP_summary.R \
  --tx_dir ~/rdeocd/txMeth \
  --panels I,J

# Optional Jul 2021 cache (32 patients, UMAP cluster counts in df_summary3):
# --validation_rdata ~/rdeocd/txMeth/ASM_validation_Jul2021_v6.RData

Panel **4J** facets compare mean DMP counts (Fisher ORs on 2×2 tables):

| Facet | Contrast |
|-------|----------|
| Allele-specific | mono-allelic (AS) vs bi-allelic DMPs |
| Clonal | all clonal vs all subclonal DMPs |
| Clonal bi-allelic only | clonal bi-allelic vs clonal mono-allelic (AS) DMPs |

Covariates (paper Fig. 4J legend; Fisher ORs on mean DMP counts):

- **Histology**: LUAD vs LUSC (adenocarcinoma vs other).
- **Tumor stage**: `1a`/`1b` vs later stages (`pTNMStage_v8` clinical update).
- **WGD** vs non-WGD (patient-level `WGD_reseg`).
- **SBS4**: high if weight &gt; 0.35, else low.
- **Number of samples sequenced**: high if &gt; 3 RRBS regions (non-LN), else low.

OR cells are **mean DMP counts per patient** in each clinical × DMP-type stratum
(then `questionr::odds.ratio` on the 2×2 mean table, matching
`Notes_ASM_plots_14Jul2021.R` ~L1136–1185). Default **`df_summary3`** in
`ASM_cohort_summary_figure4J_refresh.RData` is **one row per patient** (37 tx100
patients): CpG-level **`flag3`** counts per sample (`multip="all"`), summed across
regions. `sub_ambig_DMPs` is 0 (ambiguous class exists only in UMAP cluster labels).
Optional `--j_table summary2` uses per-sample rows. The legacy
`ASM_validation_Jul2021_v6.RData` used UMAP cluster tallies (32 patients) and
summary2 with legacy `flag` (zero subclonal in cache).

Sample-count binning: high if &gt; 3 regions sequenced (paper legend).

Third facet default is "clonal bi-allelic vs all other DMPs" (Notes ~L1149); use
`--clo_bi_contrast clo_mono` for vs clonal mono-allelic only.

## Outputs

- `outputs/Fig4A_imprinting_overview.pdf`
- `outputs/Fig4H_UMAP_<CRUKID>.pdf` (from stored `ASM_UMAP_recovered_from_pdf.RData`)
- `outputs/Fig4I_DMP_assignments_{counts,ratio}.pdf`
- `outputs/Fig4J_clinical_odds_ratios.pdf`
