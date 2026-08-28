# Aim: A number of WGBS vs WGS SNPs have a BAF difference of 1.0,
#   despite resulting in the same genotype call. This script aims to
#   investigate the cause of this discrepancy.

# Dependencies
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(fs)
options(bitmapType='cairo')
figdir = fs::dir_create('figures/240918_genotype_accuracy/')

# Example
# Chromosome 1: 591391 ; BAFDIFF == 1 in multiple samples. Sample_5 == DTB-077
ga = data.table::fread('data/pre/240820_genotyping_accuracy/genotype_scores.tsv', nrows=1000)[sample_id == 'sample_5' & bafdiff==1]

# View the RAW BAF for DTB-077 in this dataset
wgs77 = fread('/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/results/battenberg/DTB-077/DTB-077_T_normalBAF.tab')
af77 = fread('/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/results/battenberg/DTB-077/DTB-077_N_alleleFrequencies_chr1.txt')
bs77 = fread('/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/archive/2023_data/230916_wgbs_ascat/wgbs_bb_tumor_normal/results/DTB-077/Allelecounts/N/DTB-077.N.SNPs.CpGs.all.sorted.csv.gz')

# Select
wsel = wgs77[Chromosome == 1 & Position == 591391]
afsel = af77[`#CHR` == 1 & `POS` == 591391]
bsel = bs77[chrom == 1 & POS == 591391]

# Plot a subset of sites to ensure that the BAFs aren't randomised
p_wg = wgs77 %>% dplyr::sample_n(1000) %>% pull(3) %>% qplot() + labs(title='WGS BAF')
p_bs = bs77 %>% dplyr::sample_n(1000) %>% pull(BAF) %>% qplot() + labs(title='WGBS BAF')
p_baf = p_wg + p_bs
ggsave(path(figdir, 'wgbs_wgs_baf_comparison.png'), p_baf, width=20, height=10)

# Example function to return BAFs from the allelecounts file and Battenerg reference
renv::load("~/projects/wgbs_camdac_dev")
library(fs)
library(glue)

# Define allele_data before the function
get_allele_data <- function(xv) {
  # Usage: allele_data = get_allele_data('data/meta/battenberg_hg38_allele_data.csv.gz')
  allele_data=data.table::fread(xv)
  names(allele_data)=c("#CHR", "POS", "ref_sel", "alt_sel")
  allele_data$POS= as.numeric(allele_data$POS)
  setkey(allele_data, "#CHR", "POS")
  return(allele_data)
}
allele_data = get_allele_data('data/meta/battenberg_hg38_allele_data.csv.gz')

# Function to process sample data
generate_unrandomised_battenberg_baf <- function(sample_id, allele_data) {
  normal_freq = rbindlist(lapply(fs::dir_ls(glue('/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/results/battenberg/{sample_id}', sample_id=sample_id), glob='*_N_alleleFrequencies_chr*'), data.table::fread))
  setkey(normal_freq, "#CHR", "POS")
  
  norm = merge(normal_freq, allele_data, by=c("#CHR", "POS"), all.x=T)
  norm$ref_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$ref_sel)+2))])
  norm$alt_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$alt_sel)+2))])
  norm$BAF = norm$alt_count / (norm$ref_count + norm$alt_count)
  norm_table = norm[, .(Chromosome = `#CHR`, Position = POS, BAF)]
  norm_table = norm_table[!is.na(BAF)]
  names(norm_table)[3] = sample_id
  
  outdir = fs::dir_create(glue('data/out/wgs_baf/'))
  data.table::fwrite(norm_table, glue("{outdir}/{sample_id}_wgs_BAF.csv.gz"))
}

# 10 Samples for mCRCP analysis
## Note,074 is unavailable simply because the AF files were corrupted. The original is available so I genreated this manually.
## /nemo/project/proj-vanloo-secure/working/nmensah/data/210916_SU2C_WGS/allele_counts/DTB-074_N.alleleCount.txt
nkeep = readr::read_lines('~/data/210828_SU2C_normals/10_matched.txt')
# Generate the unrandomised BAF using parallele for all samples in nkeep
library(foreach)
library(doParallel)
doParallel::registerDoParallel(10)
foreach::foreach(i=nkeep, .combine='rbind') %dopar% {
    print(i)
  generate_unrandomised_battenberg_baf(i, allele_data)
}
doParallel::stopImplicitCluster()

# N.B. : Code to generate allele_data
#chr_names = c(1:22, 'X', 'Y')
#G1prefix = '/camp/home/mensahn/proj-vanloo/reference_files/human/references/Battenberg_hg38/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr'
#normal_input_data = Battenberg:::concatenateAlleleCountFiles(AFprefix, ".txt", chr_names)
#allele_data = Battenberg:::concatenateG1000SnpFiles(G1prefix, ".txt", chr_names)
#allele_data=data.table(allele_data)

