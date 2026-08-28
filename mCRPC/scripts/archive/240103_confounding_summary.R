# Create a proof-of-principle figure for summarising confounding across the cohort

# Libraries ----
source("src/fig1_confounding.R")
library(foreach)


# Inputs ----

# Pure methylation objects (.qs from 2022 CAMDAC)
fmeths <- fs::dir_ls(
    "/nemo/project/proj-vanloo-secure/working/nmensah/data/220819_PON_ANALYSIS/camdac_validation/",
    recurse = T,
    glob = "*pure.qs"
)

# Load normal methylation panels for prostate
## Normal methylation from Skvortsova et al. 2019 (Cancer Cell) : 2m CpGs
ppanel <- fread(fs::path(
    "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male",
    "/prostate.m.csv.gz"
))
## Normal methylation from Li et al. 2020 (Nature) : 28m CpGs
cpanel <- fread(fs::path(
    "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male",
    "/ch_prostate_norm.m.csv.gz"
))

# Outputs ----
outdir <- fs::dir_create("data/pre/f1_confounding")
csvdir <- fs::dir_create(fs::path(outdir, "csv"))
figdir <- fs::dir_create(fs::path(outdir, "fig"))

# Main ----

# Load a pure methylation files and run analysis
doParallel::registerDoParallel(5)
foreach(fpure = fmeths[1:5]) %dopar% {
    pid <- extract_pure_meth_fid(fpure)
    pure <- qs::qread(fpure)

    # Select DNAme sites that are unmeth in normal
    pum <- pure[m_n <= 0.05, ]
    phm <- pure[m_n >= 0.05, ]

    pipeline_confounding(pum, phm, pid, csvdir, figdir)
}
doParallel::stopImplicitCluster()

cconf <- gather_confounding_csvs(csvdir)
fwrite(cconf, fs::path(outdir, "confounding_mcrpc_raw.csv")) # Raw data from the peak finder
# Temp: commented to preserve
#fwrite(cconf, fs::path(outdir, "confounding_mcrpc.csv")) # Manual processed data for plotting

# Post manual-editing
clnconf = fread(fs::path(outdir, "confounding_mcrpc.csv"))
maxconf <- plot_conf_cohort(clnconf, peaktype = "max")
minconf <- plot_conf_cohort(clnconf, peaktype = "min")
pconf = maxconf / minconf
ggsave(fs::path(outdir, "confounding_mcrpc_max.png"), pconf, width = 10, height = 10)
