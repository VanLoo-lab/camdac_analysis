# Figure S5

Bulk and purified methylation histograms stratified by allele-specific copy
number and purity (extended TRACERx / related analyses).

Legacy script relative to `figure1B-D_S1G-H/scripts/plot_mb_with_CN_and_purity.R`
(hard-coded paths and looping over the cohort). Prefer the figure1 CLI version
for new TRACERx / mCRPC runs when possible.

## TRACERx

Requires clinical metadata and per-patient methylation / ASCAT.m outputs under
`~/rdeocd/txMeth`. Patient focus in-script defaults include `CRUK0062`; the script
then iterates TRACERx patients for multi-sample panels.

```bash
Rscript scripts/plot_mb_with_CN_and_purity.R
```

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figureS5` (`tx_dir` uses `~/rdeocd/txMeth/`).
