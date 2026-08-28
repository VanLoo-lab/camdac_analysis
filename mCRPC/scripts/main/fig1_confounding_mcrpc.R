# Create a proof-of-principle figure for summarising confounding across the cohort

# Libraries ----
source("src/fig1_confounding.R")
library(foreach)
library(matrixStats)
devtools::load_all("~/projects/VanLoo_CAMDAC")

# Inputs ----

# QC table for MCRPC samples
qc_mcrpc <- fread('data/out/qc/mcrpc_tumor.tsv')

# Choose sample for analysis, ensuring it passed QC
sample_id <- "DTB-170"
stopifnot(sample_id %in% qc_mcrpc[tPASS==T]$patient)

# Pure methylation objects from MCRPC run
fmeths <- fs::dir_ls(
    "data/source/mcrpc/main/",
    recurse = T,
    glob = "*pure.csv.gz"
)

# Load normal methylation manifest
nm_meta <- fread('data/source/manifest_wgbs_ac_normal.tsv') %>% dplyr::filter(tissue=="prostate")
nm_obs <- fs::dir_ls('data/pre/val_normal_meth', glob="*prost*")

# Load prostate normal matrix
pnorm = qs::qread(nm_obs[[1]])


# Create prostate panel by taking the median methylation and counts
panel_m = rowMedians(as.matrix(pnorm$meth), na.rm=T)
panel_c = rowMedians(as.matrix(pnorm$cov), na.rm=T)
panel_c = round(panel_c, 0)
panel_samples = rowSums(as.matrix(!is.na(pnorm$meth)))
ppanel = cbind(pnorm$pos, panel_m, panel_c, panel_samples)
setnames(ppanel, "panel_m", "m_n")
setnames(ppanel, "panel_c", "cov_n")
setnames(ppanel, "panel_samples", "samples_n")
# Add HDI annotation
pp_M = round(ppanel$m_n * ppanel$cov_n, 0)
pp_UM = ppanel$cov_n - pp_M
ppanel_hdi = calculate_counts_hdi(pp_M, pp_UM, n_cores=10)
names(ppanel_hdi) = c("m_n_low", "m_n_high")
ppanel = cbind(ppanel, ppanel_hdi)

# Load normal methylation from Li et al. 2020 (Nature) : 28m CpGs
cpanel <- fread(fs::path(
    "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male",
    "/ch_prostate_norm.m.csv.gz"
))
cpanel_hdi = calculate_counts_hdi(cpanel$M, cpanel$UM, n_cores=10)
names(cpanel_hdi) = c("m_n_low", "m_n_high")
cpanel = cbind(cpanel, cpanel_hdi)
names(cpanel) = c("chrom", "start", "end", "M", "UM", "m_n", "m_n_low", "m_n_high")
cpanel$samples_n = NA
cpanel$cov_n = cpanel$M + cpanel$UM
names(ppanel)
cpanel = cpanel[, .(chrom, start, end, m_n, cov_n, samples_n, m_n_low, m_n_high)]

# Set pure methylation object
fpure = fmeths %>% purrr::keep(~str_detect(.x, sample_id))

# Outputs ----
outdir <- fs::dir_create("data/pre/f1_confounding/mcrpc/single/f1_dtb170")
fs::file_copy(fpure, fs::path(outdir, fs::path_file(fpure)))
fwrite(nm_meta, fs::path(outdir, "prostate_panel_meta.csv"))
fwrite(ppanel, fs::path(outdir, "prostate_panel.csv.gz"))
fwrite(cpanel, fs::path(outdir, "prostate_chinese_panel.csv.gz"))

outdir_s2 <- fs::dir_create("data/pre/f1_confounding/mcrpc/single/f1_dtb097")
fpure2 = fmeths %>% purrr::keep(~str_detect(.x, "DTB-097"))
fs::file_copy(fpure2, fs::path(outdir_s2, fs::path_file(fpure2)))
# Main ----

# Define sites hyer-methylated and un-methylated in the prostate normal panel
#   This is because we lack a matched normal. Otherwise, we could just use the m_n_o (origin)
# In bulk, HDI99 < 0.2 occurs at IQR of < 0.13 and mean cov of ~40X
# Therefore here, we filter this analysis at DNAme <= 0.1 and cov >= 30X
# The DNAme values are flipped for the hyper-methylated sites
norm_um = cpanel[ M+UM >= 30 & m <= 0.1]
norm_hm = cpanel[ M+UM >= 30 & m >= 0.9]

# Load a pure methylation files and run analysis
doParallel::registerDoParallel(5)

# Select DNAme sites that are unmeth in normal
# Note I'm using only end as this the coordinates of the chinese prostate atlas
pum <- pure[norm_um[, .(chrom, end)], on = .(chrom, end)][!is.na(CN)]
phm <- pure[norm_hm[, .(chrom, end)], on = .(chrom, end)][!is.na(CN)]

# Add coverage filter to sites >10X
pum <- pum[cov >= 10]
phm <- phm[cov >= 10]

# Save data
saveRDS(
    list(
        norm_unmeth = pum,
        norm_meth = phm
    ),
    fs::path(outdir, paste0("fig1_confounding_mcrpc_", sample_id, ".rds"))
)
