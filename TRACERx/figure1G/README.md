# Figure 1G

Direct comparison of ASCAT(.m) (RRBS) and WGS allele-specific copy number
estimates. Segments are rounded to integer copy number.

Colouring (total = nA+nB, minor = nB):
- agree — minor allele pale grey, total copy number dark grey
- disagree (minor) — yellow (RRBS), baby blue (WGS)
- disagree (total) — orange (RRBS), dark blue (WGS)

Extracted from `Notes_RRBS_WGS_Genotype_LogR_BAF_comparison.R` (CN block).

Requires TRACERx RRBS ASCAT.m and matched WGS ASCAT outputs under
`~/rdeocd/txMeth/<patient_id>/`.

## TRACERx

```bash
Rscript scripts/compare_ASCAT_RRBS_vs_WGS_CN.R \
  --work_dir ~/data/CAMDAC_paper/figure1G \
  --output_dir outputs \
  --tx_dir ~/rdeocd/txMeth \
  --patient_id CRUK0062 \
  --sample_id R1,R3,R4 \
  --reseg TRUE
```

Defaults match CRUK0062 with auto-discovery of regions that have both RRBS and
WGS ASCAT outputs, so this is equivalent when those files are present:

```bash
Rscript scripts/compare_ASCAT_RRBS_vs_WGS_CN.R
```

`--patient_id` and `--sample_id` accept comma-separated lists (`--patient_id`
prefer CRUKID; Patient IDs also accepted). PDFs are written under `--output_dir` as
`<CRUKID>_<sample>_WGS_RRBS_CN_comparison.pdf`.

Patients used in the original Notes run: CRUK0069, CRUK0062, CRUK0031
(CRUK0069, CRUK0062, CRUK0031).
