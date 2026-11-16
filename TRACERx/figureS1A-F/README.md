# Supplementary Figure S1A–F

Normal and bulk-tumour methylation-fraction densities for **CRUK0062**-R3,
HDI99 demonstration, and copy-number zooms for CpGs confidently unmethylated
in the matched normal. Output filenames use CRUKID (paths accept CRUKID or Patient).

| Panel | Content |
|-------|---------|
| A | All CpGs: mn vs mb density |
| B | Simulated HDI99 at example unmethylated CpGs |
| C | Confidently unmethylated in normal (HDI99 ⊆ [0, 0.2]) |
| D | Zoom of C (gained methylation, all copy numbers) |
| E | As D, multiplicity 2+0 (expected clonal-hyper vlines) |
| F | As D, multiplicity 2+2 |

Extracted from `Notes_plot_methylation_99_HDI.R`. Panel E fixes a Notes bug
where expected-interval vlines were drawn before `nA`/`nB` were set for 2+0.

Requires `purified_tumour_DMP_calls.RData` under
`~/rdeocd/txMeth/<patient>/Methylation/<sample>/[re_segmented/]`.

## TRACERx (defaults = paper sample)

```bash
Rscript scripts/plot_S1A_F_mb_mn_densities.R \
  --work_dir ~/data/CAMDAC_paper/figureS1A-F \
  --output_dir outputs \
  --tx_dir ~/rdeocd/txMeth \
  --patient_id CRUK0062 \
  --sample_id R3 \
  --reseg TRUE \
  --panels A,B,C,D,E,F
```

(`--patient_id` may be `CRUK0062` or `CRUK0062`; files are written as
`CRUK0062_R3_S1*.pdf`.)

Equivalent when defaults apply:

```bash
Rscript scripts/plot_S1A_F_mb_mn_densities.R
```

Purity for E–F vlines is taken from `rho_reseg` in the TRACERx clinical table
when available; otherwise the Notes default `0.73` is used. Override with
`--purity`.

PDFs are written under `--output_dir` as
`<patient>_<sample>_S1{A–F}_*.pdf`.

## Related
S1G–H (bulk histograms × CN) live in `../figure1B-D_S1G-H/`.
