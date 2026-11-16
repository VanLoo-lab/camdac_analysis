# Figure 3D–E

Compare pure-tumour (mt–mt) versus bulk (mb–mb) differentially methylated
position (DMP) overlap for TRACERx, including a focused intra-tumour analysis
of CRUK0062 (CRUK0062) and genome-wide inter-/intra-patient sampling by purity.

Main runnable scripts write under `outputs/`. The exploratory PhD notebook
(`Notes_mt_mt_vs_mb_mb_DMPs_figure5FG.R`) lives under `~/data/CAMDAC_paper/figure3DE/scripts/`.

## CRUK0062 / CRUK0062 intra-tumour DMP overlap

Requires per-region purified DMP calls under `~/rdeocd/txMeth/CRUK0062`:

```bash
Rscript scripts/caculate_dmp_overlap_cruk0062.R
```

Outputs (Euler / Venn-style PDFs and related objects) are written to
`~/data/CAMDAC_paper/figure3DE/outputs`.

## Inter- and intra-tumour DMP overlap sampling

Requires genome-wide CAMDAC results
(`~/rdeocd/txMeth/CpG_wise_CAMDAC_results_noNAs.fst`) and clinical metadata:

```bash
Rscript scripts/calculate_dmp_overlap_sampling_inter_and_intra_tumour_pairs.R
```

Outputs are written to
`~/data/CAMDAC_paper/figure3DE/outputs/dmp_overlap_sampling`.
