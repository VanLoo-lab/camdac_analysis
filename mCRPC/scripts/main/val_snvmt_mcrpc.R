# Setup ----
setwd("/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac")
renv::load()
library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(fs)
library(qs)
library(foreach)
library(doParallel)
library(logger)
# Load function for reading CAMDAC-WGBS snvmt directories (WIP pipeline adadpted from WGBS-ASM)
source('src/snvmt_functions.R') 
source('src/val_snvmt.R')
options(bitmapType='cairo')
outfile = fs::path("data/pre/f4_snvmt/mcrpc_tumor.tsv")
qcfile = fs::path("data/out/qc/mcrpc_tumor.tsv")
outdir = dir_create("data/pre/f4_snvmt/results")
outdir_root = fs::path_dir(outdir)
outdir_comb = dir_create(fs::path(outdir_root, "combined"))
figdir = dir_create("data/pre/f4_snvmt/figures")
snvph= "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/snvph"
snvdp = "/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/results/snv_dpc_txt" 
main= "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/main"


# Main ----

# Get files with SNV phasing data
indirs = fs::dir_ls(snvph, recurse=T, type="d", glob="*SNV_phase/T")

# Filter by QC pass samples
qcres = fread(qcfile)[tPASS==T]$patient
indirs = indirs[match(qcres,str_extract(indirs, "DTB-\\d+"))]

# Load snvmt data for current sample
doParallel::registerDoParallel(cores=4)
res = foreach(curdir=indirs, .errorhandling="pass") %dopar% {

    log_info("Processing {curdir}")
    sid = get_sid(curdir)
    valphas = load_snvmt_data(curdir, snvdp, main)
    
    # Filter for high ccf
    valphas = valphas[hap_ccf>.95]
    valphas = valphas_to_minimal(valphas)

    # Save data
    data_outfile = fs::path(outdir, paste0(sid, "_valphas.csv.gz"))
    log_info("Saving data to {data_outfile}")
    fwrite(valphas, data_outfile)

    ## Create plot function
    log_info("Creating plots for {sid}")
    gvaph = plt_snv_phased(valphas, "hap_BAF", "snv_VAF", "VAF")
    gvmb = plt_snv_phased(valphas, "altm", "m", "Bulk")
    gvmt = plt_snv_phased(valphas, "altm", "m_t", "Pure")
    p_valphas = patchwork::wrap_plots(gvaph, gvmb, gvmt, nrow=1) & plot_annotation(title=sid)
    plt_outfile = fs::path(figdir, paste0(sid, "_valphas.png"))
    ggsave(plt_outfile, p_valphas, width=10, height=4)

    log_info("Done with {sid}")
    return(c(sid, T))
}

qs::qsave(res, fs::path(fs::path_dir(outdir), "res.qs"))

# Combine results into a single output file, one with all data (high CCF)
resfiles <- rbindlist(
    lapply(fs::dir_ls(outdir, glob="*valphas.csv.gz"), function(x){
        patient = str_extract(basename(x), "DTB-\\d+")
        x = fread(x)[, patient:=patient]
        return(x)
    })
)
fwrite(resfiles, fs::path(outdir_comb, "all_valphas.csv.gz"))

# Generate a figure for the combined data
gvaph = plt_snv_phased(resfiles, "hap_BAF", "snv_VAF", "VAF")
gvmb = plt_snv_phased(resfiles, "altm", "m", "Bulk")
gvmt = plt_snv_phased(resfiles, "altm", "m_t", "Pure")
p_valphas = patchwork::wrap_plots(gvaph, gvmb, gvmt, nrow=1) & plot_annotation(title="Combined")
plt_outfile = fs::path(figdir, "all_valphas.png")
ggsave(plt_outfile, p_valphas, width=10, height=4)

# Apply filters and save
rescln = resfiles[cov_t>=20,]
fwrite(rescln, fs::path(outdir_comb, "all_valphas_cln.csv.gz"))

# Generate a figure for the combined data
gvaph = plt_snv_phased(rescln, "hap_BAF", "snv_VAF", "VAF")
gvmb = plt_snv_phased(rescln, "altm", "m", "Bulk")
gvmt = plt_snv_phased(rescln, "altm", "m_t", "Pure")
p_valphas = patchwork::wrap_plots(gvaph, gvmb, gvmt, nrow=1) & plot_annotation(title="Combined_20Xcovt")
plt_outfile = fs::path(figdir, "all_valphas_cln.png") 
ggsave(plt_outfile, p_valphas, width=10, height=4)

# Generate a figure showing the distance from the diagonal for the methylation data
#   Here, this is simply the aboslute distance between the two values
#   This is calculated per patient from the phased snv data that we plot above.
d_nrpcc = fread(qcfile)[, .(patient, nrpcc, nrpcc_d)]
resfile_dist = rbind(calc_absdiff(resfiles, "altm", "m"), calc_absdiff(resfiles, "altm", "m_t")) %>%
    merge(d_nrpcc, by="patient")
rescln_dist = rbind(calc_absdiff(rescln, "altm", "m"), calc_absdiff(rescln, "altm", "m_t")) %>%
    merge(d_nrpcc, by="patient")
p_ndiag_resfile = plt_ndiag(resfile_dist)
ggsave(fs::path(figdir, "all_valphas_ndiag.png"), p_ndiag_resfile, width=10, height=4)
p_ndiag_rescln = plt_ndiag(rescln_dist)
ggsave(fs::path(figdir, "all_valphas_cln_ndiag.png"), p_ndiag_rescln, width=10, height=4)

# Additionall, compare the empirical nrpcc with the estimate 
plt_nrpcc = ggplot(resfile_dist, aes(x=e_nrpcc, y=nrpcc)) +
    geom_point() +
    geom_abline(intercept=0, slope=1) +
    theme_classic() +
    labs(x="phased_snv_nrpcc")
ggsave(fs::path(figdir, "all_valphas_nrpcc.png"), plt_nrpcc, width=4, height=4)
