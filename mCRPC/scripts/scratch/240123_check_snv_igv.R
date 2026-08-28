# Check SNVs in IGV that have poor overlap between WGBS VAF and WGS VAF
library(data.table)
library(logger)
library(fs)
source('src/val_snvmt.R')
options(bitmapType='cairo')

# Meta file if needed
srr = fread('data/source/table_s1_mcrpc_srr.csv')
srr[submitted_subject_id=="DTB-034"][]

# Load snv files
snvph = fs::path('/camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/pre/f4_snvmt/results')
snvphf = dir_ls(snvph)


# Focusing on DTB-034
# /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/pre/f4_snvmt/figures/DTB-034_valphas.png
VF = fread(snvphf %>% purrr::keep(~str_detect(.x, "DTB-034")))
VF[, diff:=hap_BAF - snv_VAF]
boxplot.stats(VF$diff)$stats
VF[snv_REF == "C" & snv_ALT == "G", ][chrom==5 & hap_POS == 58736451]

# Opening IGV to review
## Input WGS : open smb://mensahn@data2.thecrick.org/proj-vanloo-secure/working/nmensah/data/210916_SU2C_WGS/bam
## Input WGBS : open smb://mensahn@data2.thecrick.org/proj-vanloo/secure/working/nmensah/su2c_ac/bam_symlink
VF2 = VF[, .(md=mean(diff)), by=c("snv_REF", "snv_ALT")] %>% dplyr::arrange(desc(md))


## Test another one just for sanitty
XF2 = fread(snvphf[37]) %>% dplyr::mutate(diff=hap_BAF - snv_VAF) %>% dplyr::group_by(snv_REF, snv_ALT) %>% 
    dplyr::summarise(md=mean(diff)) %>% dplyr::arrange(desc(md)) %>% data.table
XF2

## Plot without the failers
resfiles = fread('data/pre/f4_snvmt/combined/all_valphas.csv.gz')
rf2 = resfiles[
    !(snv_REF=="G" & snv_ALT=="C"),][
    !(snv_REF=="T" & snv_ALT=="G"),][
    !(snv_REF=="C" & snv_ALT=="G"),][
    !(snv_REF=="A" & snv_ALT=="C"),]
XF3 = rf2 %>% dplyr::mutate(diff=hap_BAF - snv_VAF) %>% dplyr::group_by(snv_REF, snv_ALT) %>% 
    dplyr::summarise(md=mean(diff)) %>% dplyr::arrange(desc(md)) %>% data.table

gvaph = plt_snv_phased(rf2, "hap_BAF", "snv_VAF", "VAF")
figdir = "figures/241018_snv_phase_check"
ggsave(fs::path(figdir, "errSNV_all_valphas.png"), gvaph, width=5, height=5)

gvaph = plt_snv_phased(rf2[cov_t>=20], "hap_BAF", "snv_VAF", "VAF")
figdir = "figures/241018_snv_phase_check"
ggsave(fs::path(figdir, "errSNV_all_valphas_20x.png"), gvaph, width=5, height=5)
