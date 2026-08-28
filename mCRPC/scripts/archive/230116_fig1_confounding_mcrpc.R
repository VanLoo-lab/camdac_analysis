# Create a proof-of-principle figure for summarising confounding across the cohort

# Libraries ----
source("src/fig1_confounding.R")
library(foreach)

# Inputs ----

# QC table for MCRPC samples
qc_mcrpc <- fread('data/out/qc/mcrpc_tumor.tsv')

# Pure methylation objects from MCRPC run
fmeths <- fs::dir_ls(
    "data/source/mcrpc/main/",
    recurse = T,
    glob = "*pure.csv.gz"
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
outdir <- fs::dir_create("data/pre/f1_confounding/mcrpc")
csvdir <- fs::dir_create(fs::path(outdir, "csv"))
figdir <- fs::dir_create(fs::path(outdir, "fig"))
ctable_qc <- fs::path(outdir, "confounding_mcrpc.csv") # confounding table for manual QC edits
conf_dot_plot <- fs::path(outdir, "confounding_mcrpc_dot.png")
conf_tile_plot <- fs::path(outdir, "confounding_mcrpc_tile.png")

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
foreach(fpure = fmeths) %dopar% {

    # Get patient identifier
    pid <- extract_pure_meth_fid(fpure)
    print(sprintf("Processing %s", pid))

    # Skip if data exists
    csvhm <- fs::path(csvdir, paste0(pid, "_HM.csv"))
    csvnu <- fs::path(csvdir, paste0(pid, "_NU.csv"))

    if (all(fs::file_exists(c(csvhm, csvnu)))) {
        log_info("CSVs already exist. Skipping.")
        return()
    }


    # Load data
    pure <- data.table::fread(fpure, select=c("chrom", "start", "end", "CN", "p", "m", "cov"))

    # Select DNAme sites that are unmeth in normal
    # Note I'm using only end as this the coordinates of the chinese prostate atlas
    pum <- pure[norm_um[, .(chrom, end)], on = .(chrom, end)][!is.na(CN)]
    phm <- pure[norm_hm[, .(chrom, end)], on = .(chrom, end)][!is.na(CN)]

    # Add coverage filter to sites >10X
    pum <- pum[cov >= 10]
    phm <- phm[cov >= 10]

    pipeline_confounding(pum, phm, pid, csvdir, figdir)
}
doParallel::stopImplicitCluster()

stop() # Stop here to QC the plots

# Gather CSVs, filtering out any records where QC==1 (failed to find peak)
cconf <- gather_confounding_csvs(csvdir)

# Filter to samples that passed QC
qcpass = qc_mcrpc[tPASS==T,]$patient
# Note: not all samples will have peaks found (n=82/92 with pure meth)
cconf_id = str_extract(cconf$id, "DTB-\\d+") 
# Note: not all samples pass qc (n=83)
cconf = cconf[cconf_id %in% qcpass,]
# Remainder for this analysis: (n=77) 


# Note: Manually edited to use 5+ rather than 6+ as max CN state
# This is in the code now but mcrpc was run and qc initially without
cconf$CN_state = ifelse(cconf$CN_state %in% c("5", "6+"), "5+", cconf$CN_state)

# Write result if file does not exist
if (!fs::file_exists(ctable_qc)) {
    fwrite(cconf, ctable_qc)
} else {
    print(sprintf("QC table already exists at %s", ctable_qc))
}

stop()

# Generate plots ----

# Load data
clnconf <- fread(ctable_qc)

# Create dot plot
conf_dot <- plot_conf_dot(clnconf, conf_dot_plot)

# Create tile plot
conf_tile <- plot_conf_tile(clnconf, conf_tile_plot)
