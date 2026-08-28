# Aim: Count the unmethylated CpG sites in the pure methylation data for each mCRPC sampel
#  N.B. Count sites enriched for the normal methylation panel (SU2C/skvorstova) and the chinese atlas panel


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
conf_count = foreach(fpure = fmeths, .combine='rbind') %dopar% {
    pid <- extract_pure_meth_fid(fpure)
    pure <- qs::qread(fpure)

    # Select DNAme sites that are unmeth in normal
    pum <- nrow(pure[m_n<=0.05 & m <= 0.05,]) / nrow(pure[m_n<=0.05,])   
    phm <- nrow(pure[m_n>=0.05 & m >= 0.05,]) / nrow(pure[m_n>=0.05,])   
    ahm <- nrow(pure)
    dataset <- data.table(
        pid = pid,
        numeth = pum,
        nhyper = phm,
        allsites = ahm
    )
    
    return(dataset)

}
doParallel::stopImplicitCluster()

outdir = fs::dir_create('data/pre/f1_confounding/')
fwrite(conf_count, fs::path(outdir, "confound_count.csv")) 

# Summary stats
round(mean(conf_count$numeth),2)

