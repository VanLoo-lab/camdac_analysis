# Figure 1H–K

Compare genomic (WES/WGS) and methylation-derived (RRBS/WGBS) allele-specific
copy number, purity, and ploidy across TRACERx, mCRPC, and CRC.

`rrbs_dnaseq_intersections.R` builds TRACERx RRBS–DNA-seq ASCN overlap inputs.
`figure1H-K.R` loads those (and other cohort overlap tables under `inputs/`)
and writes the figure panels to `outputs/`.

## TRACERx ASCN intersections (upstream)

Requires access to TRACERx methylation and WES ASCAT outputs:

```bash
Rscript scripts/rrbs_dnaseq_intersections.R \
  --wdir ~/data/CAMDAC_paper/figure1H-K \
  --txmeth_dir ~/rdeocd/txMeth \
  --txwes_dir /path/to/TRACERx_WES
```

TRACERx RRBS–WGS ASCN overlap for the samples with WGS ASCAT solutions
(currently n = 7) is built from RRBS + WGS outputs under `txMeth`:

```bash
Rscript scripts/rrbs_wgs_intersections.R \
  --wdir ~/data/CAMDAC_paper/figure1H-K \
  --txmeth_dir ~/rdeocd/txMeth
```

This writes `inputs/TRACERx_RRBS_WGS_overlap.txt` and
`inputs/dt_tracerx_rrbs_wgs.RData`.

## Figure panels

When run with `Rscript`, the working directory defaults to
`~/data/CAMDAC_paper/figure1H-K` (clinical/ASCAT.m paths use `~/rdeocd/txMeth`).
Precomputed overlap tables must be present under `inputs/`.

```bash
Rscript scripts/figure1H-K.R
```

Outputs are written under `outputs/`, including purity and ploidy comparisons,
histology-stratified summaries, LOH vs ploidy, LUAD/LUSC gain–loss profiles,
and BS-seq vs DNA-seq ASCN overlap.
