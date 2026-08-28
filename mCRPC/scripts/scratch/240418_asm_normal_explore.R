# Aim: Investigate ASM results at CpGs in normal samples

# Setup ----
library(data.table)
library(gridExtra)
library(logger)
library(dplyr)
library(patchwork)
options(bitmapType='cairo')
source('src/asm_icr.R')

# Inputs ----
figdir = fs::dir_create('figures/240418_asm_icr')
outdir = fs::dir_create('data/archive/2024_data/240418_asm_')
icr_f = '/nemo/project/proj-vanloo-secure/working/nmensah/data/240328_CRC_LESKY/240328_lesky_asm/240328_icrs_imprintome_plus_geneimprint.csv'
norm03 = '/nemo/project/proj-vanloo-secure/working/nmensah/data/230919_SU2C_ASCAT_MATCHED/results/DTB-077/Allelecounts/N/DTB-077.N.SNPs.CpGs.all.sorted.csv.gz'
norm24 = '/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_SU2C_CAMDAC/results/240418/DTB-077/Allelecounts/N/DTB-077.N.SNPs.CpGs.all.sorted.csv.gz'
phase_dir = '/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_SU2C_CAMDAC/results/240418/DTB-077/AlleleSpecific/N'

# Data ----
log_info("Setup")
icr = fread(icr_f) %>% dplyr::mutate(chrom=gsub('chr', '', chrom)) %>% dplyr::filter(gene_imprint_db==TRUE); setkey(icr, chrom, start, end)
dnorm03 = parse_norm(norm03) %>% foverlaps(icr, type="any", mult="all", nomatch=NA) %>% dplyr::filter(!is.na(id))
dnorm24 = parse_norm(norm24) %>% foverlaps(icr, type="any", mult="all", nomatch=NA) %>% dplyr::filter(!is.na(id))

# Figures ----
log_info("Plotting initial ICR data")
create_icr_boxplot_count(dnorm03, 'icr_normal_03.png')
create_icr_boxplot_count(dnorm24, 'icr_normal_24.png')
create_icr_median_hist(dnorm03, 'icr_normal_03_hist.png')
create_icr_median_hist(dnorm24, 'icr_normal_24_hist.png')
dicr03 = get_icr_med_window(dnorm03)
dicr24 = get_icr_med_window(dnorm24)

# Load ASM data and filter to ICR regions with ASM signal in the normal
log_info("Loading ASM data")
asm = helper_load_asm(phase_dir, chrom="all", cln_output=T)
asm_filt = unique(asm %>% dplyr::select(-qname)) %>% dplyr::filter(!is.na(altm) & !is.na(refm))
asm_cln = foverlaps(asm_filt, icr[id %in% unique(dicr24$id)], type="any", mult="all", nomatch=NA) %>% dplyr::filter(!is.na(id))

# Plot difference in ASM between alleles per ICR 
log_info("Plottig ASM data")
p_icr = create_icr_plot(asm_cln) + labs(title="DTB-077-N-ASM")
ggsave(fs::path(figdir, 'icr_normal_asm_diff.png'), p_icr, width=10, height=5)

# Summary stat: for each ICR, tally the number of CpGs and the percentage of those
#   with an abs(alt-ref) meth >= 0.25 and >= 0.5
log_info("Plotting ASM stats")
asm_stat = calculate_asm_statistics(asm_cln)
g = tableGrob(asm_stat)
ggsave(fs::path(figdir, 'icr_normal_asm_stat.png'), g, width=10, height=5)
