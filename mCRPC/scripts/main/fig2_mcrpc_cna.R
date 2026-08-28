# Script to compare CNA solutions from CAMDAC WGBS and WGS for mCRPC samples
library(data.table)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
options(bitmapType='cairo')

# Functions -----
source('src/fig2_cna.R')

# Inputs ----
# Output directory for figures
figdir = fs::dir_create(fs::path("figures/pre/f2_mcrpc_cna"))
# Output directory for results
resdir = fs::dir_create(fs::path("data/pre/f2_mcrpc_cna"))
# Root directory for CAMDAC files. Searches here for Battenberg files: _subclones.txt, _purity.txt, _ploidy.txt
wgbs_bb_dir = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/mcrpc/tusnp"
# Root directory for Battenberg. Searches here for Battenberg files: _subclones.txt, _purity.txt, _ploidy.txt 
wgs_bb_dir = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/mcrpc/wgs_battenberg"


# Generate CNA data for both cohorts ----
wgbs_segments = fs::path(resdir, "wgbs_bb_tusnp_segments.csv.gz")
aggregate_battenberg_cna_segments(wgbs_bb_dir, outfile=wgbs_segments)
wgs_segments = fs::path(resdir, "wgs_bb_segments.csv.gz")
aggregate_battenberg_cna_segments(wgs_bb_dir, outfile=wgs_segments)

# Optional:: Filter to samples that passed QC ----
qc_pass = fread("data/out/qc/mcrpc_tumor.tsv")[tPASS==T]$patient
wgbs_segments_qc = fread(wgbs_segments)[ grepl(paste0(qc_pass, collapse="|"), sample_id) ]
wgs_segments_qc = fread(wgs_segments)[ grepl(paste0(qc_pass, collapse="|"), sample_id) ]
# Overwrite segments file
fwrite(wgbs_segments_qc, wgbs_segments)
fwrite(wgs_segments_qc, wgs_segments)

# Get WGS and WGBS intersect data for a cohort using segments files
id = get_intersect_data(wgs_segments, wgbs_segments, "WGS", "WGBS")

# Write result
fwrite(id$res, fs::path(resdir, "mcrpc_wgbs_wgs_overlap.csv"))
qs::qsave(id$intersects, fs::path(resdir, "mcrpc_wgbs_wgs_interects.qs"))

# Plot intersect data
plot_intersect_data(id$res, fs::path(figdir, "mcrpc_cna_validation.png"))

print("Done")
