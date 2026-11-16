# Figure 3A–C (also S7D–F)

DMP calling performance on **simulated** data stratified by copy number and
purity.

| Panel | Content |
|-------|---------|
| 3A | Bi-allelic tumor–normal **FNR** vs copy number |
| 3B | Bi-allelic tumor–normal **FPR × 10³** vs copy number |
| 3C | Bi-allelic tumor–tumor **FNR/FPR** CN heatmaps |
| S7D | Mono/bi tumor–normal **FNR** vs copy number |
| S7E | Mono/bi tumor–normal **FNR** vs coverage interval |
| S7F | Mono-allelic tumor–tumor **FNR/FPR** CN heatmaps |

**Source:** `R_scripts/Notes_DMP_simulation.R`.

### Regenerate Figure 3A–C and S7D–F

```bash
Rscript scripts/plot_3A_C_DMP_simulation.R \
  --work_dir ~/data/CAMDAC_paper/figure3A-C \
  --sim_rdata ~/rdeocd/txMeth/power_sim_results_list.RData \
  --panels A,B,C,D,E,F
```

The tumor–normal beta comparisons are accelerated with Rcpp. Figure 3A–C is
fixed to the bi-allelic class; S7D/E compare mono- and bi-allelic simulations,
and S7F is the mono-allelic counterpart to 3C.

Main outputs:

- `Fig3A_TN_FNR_vs_copy_number.pdf`
- `Fig3B_TN_FPR_vs_copy_number.pdf`
- `Fig3C_TT_FNR_heatmap_2x3.pdf` (pale blue → dark orange)
- `Fig3C_TT_FPR_heatmap_2x3.pdf` (blue → red)
- `Fig3C_TT_FNR_FPR_heatmap_2x6.pdf` (panel C layout)
- `Fig3A-C_simulated_DMP.pdf` (assembled layout)
- `FigS7D_TN_FNR_vs_copy_number.pdf`
- `FigS7E_TN_FNR_vs_coverage_interval.pdf`
- `FigS7F_TT_mono_FNR_FPR_heatmap_2x6.pdf`
- `FigS7D-F_simulated_DMP.pdf` (assembled supplementary layout)
- source-rate TSVs for each panel

S7F plots only the mono-allelic CN states 2, 4, and the capped `>=6` bin.

Filenames use CRUKID only.
