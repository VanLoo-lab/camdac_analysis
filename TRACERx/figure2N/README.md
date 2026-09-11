# Figure 2N

Pairwise Pearson correlation distances between bulk, CAMDAC-purified, and
FACS-purified methylation profiles across TRACERx samples.

Builds / reuses matrices under `inputs/` (`lung_mb.RDS`, `lung_mt.RDS`,
`lung_mf.RDS`) and writes distance plots to `outputs/`.

## TRACERx

Requires TRACERx RRBS (`~/rdeocd/txMeth`) and FACS-sorted methylation data
(`~/rdeocd/FACS_sorted_txMeth`):

```bash
Rscript scripts/plot_pearson_distances.R
```

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figure2N`.
