# Figure 3F–G

Evaluate CAMDAC DMR recall at known LUAD / LUSC promoter DMRs (Belinsky gene
set) and compare with bulk tumour calls.

Uses promoter lists and precomputed DMR tables under `inputs/`, plus TRACERx
clinical metadata under `~/rdeocd/txMeth`. Writes recall plots and tables to
`outputs/`.

## TRACERx

```bash
Rscript scripts/compute_recall_of_known_dmr.R
```

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figure3FG` (`tx_dir` uses `~/rdeocd/txMeth`).
