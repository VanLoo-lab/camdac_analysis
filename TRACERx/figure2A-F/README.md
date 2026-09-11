# Figure 2A–F

Compare naive (bulk) and CAMDAC methylation rates for TRACERx normals and
tumours, including SNP-context stratifications.

Requires TRACERx clinical metadata and per-patient allele-count / methylation
objects under `~/rdeocd/txMeth`. Outputs are written under `outputs/`.

## TRACERx

```bash
Rscript scripts/Comparing_naive_and_CAMDAC_methylation.R
```

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figure2A-F` (`tx_dir` uses `~/rdeocd/txMeth/`).
