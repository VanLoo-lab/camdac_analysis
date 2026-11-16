# Supplementary Figure S3A–D

RRBS heterozygous SNP genotyping accuracy vs WGS ground truth (TRACERx pilot
normals). SNPs are heterozygous when `0.1 ≤ BAF ≤ 0.9`.

| Panel | Content |
|-------|---------|
| A | FP rates by SNP type / sequence context (top) + BAF-error histograms at WGS-homozygous SNPs (bottom) |
| B | FN rates by SNP type / sequence context (top) + BAF-error histograms at WGS-heterozygous SNPs (bottom) |
| C | FNR vs FPR across minimum coverage and heterozygosity boundaries |
| D | TP/FN enrichment at near-het-CCGG / short-read / low mapping-quality loci (**CRUK0062**) |

**Source:** `R_scripts/Notes_RRBS_WGS_Genotype_LogR_BAF_comparison.R`
(block `#paper figures 3 and S3`).

Plots from the cached pilot object
`~/rdeocd/txMeth/all_pilot_patients_genotype_TNR_FNR.RData`
(patients **CRUK0069**, **CRUK0062**, **CRUK0031**). Does **not** re-run the WGS comparison.

S3E–F (WGBS / CRC) are under
`camdac_analysis/CRC/Figure S3 Genotyping Accuracy/` and are not this script.

## Run

```bash
Rscript scripts/plot_S3_genotype_accuracy.R \
  --work_dir ~/data/CAMDAC_paper/figureS3 \
  --genotype_rdata ~/rdeocd/txMeth/all_pilot_patients_genotype_TNR_FNR.RData \
  --panels A,B,C,D
```

### Outputs (`outputs/`)

| File | Panel |
|------|-------|
| `S3A_FPR_by_SNP_type.pdf` (+ `.tsv`) | A bubbles |
| `S3A_BAF_error_hist_WGS_hom.pdf` | A histogram |
| `S3B_FNR_by_SNP_type.pdf` (+ `.tsv`) | B bubbles |
| `S3B_BAF_error_hist_WGS_het.pdf` | B histogram |
| `S3C_FNR_vs_FPR_thresholds.pdf` (+ `.tsv`) | C |
| `S3D_TP_FN_sources_of_bias.pdf` | D (CRUK0062) |
