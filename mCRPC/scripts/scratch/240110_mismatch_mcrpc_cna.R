# A number of mcrpc files did not fit the CNA well when using tumor-only SNPs
# Here, we count how many of these required refitting.
library(data.table)
library(dplyr)

# Load CNA validation results for samples with poor overlap
setwd("~/projects/paper_camdac")
dt = fread('data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_overlap.csv')
ll = dt[perc_overlap_1 <= 75] # 15 samples with <50% overlap
idents = ll$sample_id
dti = dt %>% filter(sample_id %in% idents)

# Load mcrpc cna REFIT data
qcfile = fread('/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/files/qc_battenberg_wgs.csv')
refitfile = fread('/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/files/240115_refit_battenberg_wgs.tsv', sep='\t', header=F)
qc_refits = qcfile[refit!="",]$dtb_id
rf_refits = unique(refitfile[[1]])
frac_refit = sum(rf_refits %in% qc_refits)/length(rf_refits) # Refits
frac_all = sum(qc_refits %in% rf_refits)/length(qc_refits) # All 
# 1 Extra sample in rf_refits, likely not updated in qc_refits. Use rf as the marker
n_refit = length(rf_refits)
n_wrong = length(idents)
# Only 57% (11/18) such samples were in rf_refits (n=19)
fract_wrong_refit = sum(idents %in% rf_refits)/length(rf_refits)

# Solution: Create a slide deck where we review each CNA solution plot for the 15 samples.
# Easiest would be to give a directory, and a sample ID and create a new composite PNG from the required battenberg files.
# Root directory for CAMDAC files. Searches here for Battenberg files: _subclones.txt, _purity.txt, _ploidy.txt
wgbs_bb_dir = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/mcrpc/tusnp"
# Root directory for Battenberg. Searches here for Battenberg files: _subclones.txt, _purity.txt, _ploidy.txt 
wgs_bb_dir = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/mcrpc/wgs_battenberg"

# WARN: NO ASPCF plot for wgs

# file globs:
# *Profile_average.png
# *_distance.png
# *.tumor.p[ng]
library(grid)
library(gridExtra)
library(png)

load_bb_imgs <- function(pid, bb_dir, outfile, topt){
    # Find images
    bsregex = paste0(pid, c(".*Profile_average.png",".*T_distance.png",".*tumour.png"), collapse="|")
    bsfiles = fs::dir_ls(bb_dir, recurse=T, regex=bsregex)
    bsfiles = bsfiles[order(bsfiles)]

    # Read images
    bsgrobs = lapply(bsfiles, function(x) rasterGrob(as.raster(readPNG(x)), interpolate=F))
    layoutmat = rbind( c(1, 1, 1, 2, 2), c(3,3,3,NA,NA))
    bsplot = marrangeGrob(bsgrobs, ncol=2, nrow=2, layout_matrix=layoutmat, top=topt)
    ggsave(outfile, bsplot, width=25, height=10)
}

figdir = fs::dir_create("figures/240110_mismatch_mcrpc_cna")

library(foreach)
library(doParallel)
registerDoParallel(cores=8)
foreach(pid = idents) %dopar% {
    print(sprintf("Processing %s", pid))
    wgs_outfile = fs::path(figdir, paste0(pid, "_wgs.png"))
    topname = ifelse(pid %in% rf_refits, "WGS_REFIT", "WGS")
    load_bb_imgs(pid, wgs_bb_dir, wgs_outfile, topname)

    wgbs_outfile = fs::path(figdir, paste0(pid, "_wgbs.png"))
    load_bb_imgs(pid, wgbs_bb_dir, wgbs_outfile, "WGBS")
}

# Final comparison result:
##  023 - WGBS ploidy incorrect. Use WGS.
##  034 - WGBS ploidy incorrect. Use WGS.
##  063 - Highly aberrant. No clear solution. QC_DROP
##  069 - REFIT WGS (not WGD). Use WGS. [In Progress] ::
##  074 - REFIT WGS. Was not initially analysed. [In Progress] ::
##  089 - WGBS failed. QC_DROP
##  104 - WGBS ploidy incorrect. Use WGS.
##  111 - WGBS ploidy incorrect. Use WGS.
##  119 - WGBS ploidy incorrect. Use WGS.
##  121 - WGBS ploidy incorrect. Use WGS.
##  128 - WGBS ploidy incorrect. Use WGS.
##  141 - WGBS ploidy incorrect. Use WGS.
##  151 - WGBS ploidy incorrect. Use WGS.
##  202 - WGBS ploidy incorrect. Use WGS.
##  204 - WGBS ploidy incorrect. Use WGS.
##  205 - WGBS failed. QC_DROP
##  216 - WGBS ploidy incorrect. Use WGS.
##  223 - WGBS ploidy incorrect. Use WGS.
##  266 - REFIT WGS (WGD). plodiy 3.4, purity 0.6


# Cont> 24 Jan
# Get data with mismatched ploidies
dtmp = dt[abs(ploidy_a - ploidy_b) > .5]
dtmp
mpfigdir = fs::dir_create("figures/240124_mismatch_mcrpc_cna")

rfpid = c("DTB-104", "DTB-128", "DTB-151", "DTB-202", "DTB-204", "DTB-223")
mrfdir = fs::dir_create("figures/240125_mismatch_mcrpc_cna")

library(foreach)
library(doParallel)
library(ggplot2)
options(bitmapType='cairo')
registerDoParallel(cores=8)
foreach(pid = rfpid) %do% {
    print(sprintf("Processing %s", pid))
    wgs_outfile = fs::path(mrfdir, paste0(pid, "_wgs.png"))
    load_bb_imgs(pid, wgs_bb_dir, wgs_outfile, "WGS")

    wgbs_outfile = fs::path(mrfdir, paste0(pid, "_wgbs.png"))
    load_bb_imgs(pid, wgbs_bb_dir, wgbs_outfile, "WGBS")
}
stopImplicitCluster()

# Progress 24 Jan
# Refit most files: /nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/data/in/refit_tusnp.txt
# However, DTB-128 will likely double again on WGS: .68 and 7

# Check @1903pm

# Note: refit involves deleting the CNA directory then re-running tusnp with
# the appropriate sid
# TODO: Refits are complete but can't get the same ploidy for many samples
##  I've been informed by Maxime that there is a bb function we can override to enforce the purity/ploidy window
##  Implementing override now.
## DTB-023 OK
## DTB-034 OK
## DTB-074 OK
## DTB-104 No valid CNA solution
## DTB-111 OK
## DTB-121 OK
## DTB-128 Refitting WGBS ???
## DTB-141 OK
## DTB-151 Refitting WGBS ???
## DTB-202 Refitting WGBS ???
## DTB-204 Refitting WGBS ???
## DTB-216 OK
## DTB-223 Refitting WGBS ???

