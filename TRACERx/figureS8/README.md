# Figure S8 — Imprinted loci (TRACERx)

| Panel | Content |
|-------|---------|
| A | IGF2/H19 locus (CRUK0073-N / CRUK0073-R2; CRUK0073) |
| B | MEST/MESTIT1 locus (YORC-25-N / YORC-25-T) |
| C | NSCLC RRBS overview of LOI / LOH / Imprinted status |
| D | CRC WGBS overview of LOI / LOH / Imprinted status |

**Source:** `R_scripts/Notes_imprinted_genes_v2.R`.

### Regenerate panel C

```bash
Rscript scripts/plot_S8_imprinting.R \
  --work_dir ~/data/CAMDAC_paper/figureS8 \
  --tx_dir ~/rdeocd/txMeth \
  --panels C
```

Optional: `--imprint_rdata` to override the Oct2021 gene-subset cache;
`--panels C,C_H19` also writes an H19 sample strip.

Writes `S8C_NSCLC_imprinting_overview.pdf` (+ TSVs) with CRUK sample ids.

### Seeded panel A

Copied Notes allele-specific H19 tracks for CRUK0073 (= CRUK0073):

- `S8A_CRUK0073_N_IGF2_H19.pdf`
- `S8A_CRUK0073_R2_IGF2_H19.pdf`

Panels B (YORC-25 MEST) and D (CRC WGBS overview) need CRC inputs not
present under `txMeth/`; leave deferred until those caches are located.
