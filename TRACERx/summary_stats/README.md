# CAMDAC summary statistics

The CLI reports per-sample phased SNP/CpG counts and coverage statistics.

```bash
# Bulk bi-allelic methylation counts (default)
Rscript scripts/number_of_CpG_phased_SNPs.R \
  --tx_dir ~/rdeocd/txMeth \
  --work_dir ~/data/CAMDAC_paper/summary_stats \
  --source allelecounts

# Phased ASM allele-specific counts (writes *_ASM outputs)
Rscript scripts/number_of_CpG_phased_SNPs.R \
  --tx_dir ~/rdeocd/txMeth \
  --work_dir ~/data/CAMDAC_paper/summary_stats \
  --source asm
```

Optional filters:

```bash
--patients CRUK0034,CRUK0090 --samples N,R1,R2
```

Use `--schema_only` without accessing the data mount to print the required
allele-count fields.

## Coverage definitions

- `mean_cpg_coverage_all`: mean `total_counts_m` across unique CpGs.
- `mean_cpg_coverage_non_snp`: same, where `POS` is `NA`.
- `mean_cpg_coverage_snp`: same, where `POS` is present.
- **`--source asm` only** (from `df_phased_methylation`, not ref/alt labels):
  `mean_cpg_coverage_genotype_1` / `_2` and `_snp` / `_non_snp` variants use
  `total_counts_m_1` / `total_counts_m_2` on tumor (or `total_counts_m_n_1` /
  `total_counts_m_n_2` on normal), i.e. read counts phased to `genotype_1` and
  `genotype_2`.
- `mean_snp_depth`: mean `total_depth` across unique SNPs. This includes
  `other_counts`.
- `mean_informative_snp_depth`: mean `total_counts` across unique SNPs, where
  `total_counts = ref_counts + alt_counts`.

Counts are first collapsed by genomic locus so loci represented by multiple
records are not overweighted.

- **`--source allelecounts`**: `{Patient}/Allelecounts/{sample}/{Patient}.{sample}.SNPs.CpGs.all.sorted.RData`
  (loads `gr_normal` for `N`, `gr_tumour` for `R*`, converted to a table — not
  `Methylation/gr_*_m.RData`).
- **`--source asm`**: phased ASM methylation rates plus `allele_specific_allele_counts_*.RData`.

The `patient` column uses **CRUKID** when the clinical table is present under
`tx_dir`; paths in `allelecounts_file` still use internal Patient folders (looked up from CRUKID).

## Outputs

Bulk (no suffix):

- `outputs/unique_snps_cpgs_summary.{tsv,csv}`
- `outputs/cohort_summary.tsv`
- `outputs/unique_snps_cpgs_analysis.RData`
- `outputs/log.out`

ASM (`_ASM` suffix):

- `outputs/unique_snps_cpgs_summary_ASM.{tsv,csv}`
- `outputs/cohort_summary_ASM.tsv`
- `outputs/unique_snps_cpgs_analysis_ASM.RData`
- `outputs/log_ASM.out`
