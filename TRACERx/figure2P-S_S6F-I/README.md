# Figure 2P–S / S6F–I

Validate CAMDAC deconvolution using SNVs in tumour LOH regions, and compare
SNV-phased tumour / wild-type methylation to CAMDAC purified and adjacent-normal
estimates (TRACERx RRBS; CRC RRBS and WGBS).

`Phase_methylation_at_SNVs.R` phases methylation at SNVs from BAM / ASCAT.m
inputs (upstream; writes under per-patient TRACERx trees).
`SNV_deconvolution_validation.R` loads precomputed phased objects from
`inputs/` and writes figure panels under dated folders in `outputs/`.

## Phase methylation at SNVs (upstream)

Requires TRACERx clinical metadata, SNV calls, BAM paths, and CAMDAC outputs
under `~/rdeocd/txMeth` (paths are hardcoded in the script):

```bash
Rscript scripts/Phase_methylation_at_SNVs.R
```

## Deconvolution validation panels

Precomputed phased RDS/RData files must be present under `inputs/`:

```bash
Rscript scripts/SNV_deconvolution_validation.R
```

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figure2P-S_S6F-I`. Plots are written to
`outputs/<YYYY-MM-DD>/`.
