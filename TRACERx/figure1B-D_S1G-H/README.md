# Figure 1B–C

Bulk methylation histograms stratified by allele-specific copy number and purity.

Run once per project. Shared plotting logic is used for both; data loading and
figure parameters (`hdi_width`, `bin_width`, multiplicity selection) stay
project-specific so outputs match the previous hard-coded behaviour.

## TRACERx

```bash
Rscript scripts/plot_mb_with_CN_and_purity.R \
  --project tracerx \
  --work_dir ~/data/CAMDAC_paper/figure1B-D_S1G-H \
  --output_dir outputs \
  --patient_id CRUK0062 \
  --tx_dir ~/rdeocd/txMeth \
  --reseg TRUE \
  --clinical_info tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt
```

Defaults match the command above, so this is equivalent:

```bash
Rscript scripts/plot_mb_with_CN_and_purity.R --project tracerx
```

## mCRPC

```bash
Rscript scripts/plot_mb_with_CN_and_purity.R \
  --project mcrpc \
  --work_dir ~/data/CAMDAC_paper/figure1B-D_S1G-H \
  --output_dir outputs \
  --patient_id DTB-097,DTB-170 \
  --input_dir inputs \
  --panel prostate_panel.csv.gz
```

`--patient_id` accepts a comma-separated list. Outputs are written under `--output_dir`.
