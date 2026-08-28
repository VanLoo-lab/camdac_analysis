# A script to estimate genotyping accuracy in CAMDAC files
library(fs)
library(data.table)
library(stringr)
library(glue)
library(dplyr)

# Functions ----

# Convert allele counts to Battenberg *normalBAF table for WGS comparison. 
# Required as the default CAMDAC pipeline randomises SNPs.
allele_counts_to_battenberg_baf_tab <- function(ac_file, outfile){
    # ac_file: A CAMDAC allele counts file
    # outfile: A file to write the Battenberg NORMAL BAF table to
    x = data.table::fread(ac_file, select=c('chrom', 'POS', 'BAF', 'total_counts'))
    x = x[total_counts >= 10] # At least 10 reads in the normal
    x = x[!is.na(POS) & !is.na(BAF)]
    x = x[!duplicated(x[, .(chrom, POS)])]

    # Locate fields with SNPs and return BAF values. Third column is the file name
    res = x[, .(Chromosome=chrom, Position=POS, BAF=BAF)]
    ac_basename = fs::path_file(ac_file)
    names(res)[3] = ac_basename

    # Write out
    data.table::fwrite(res, outfile, sep="\t")
    print(paste0("Wrote ", outfile))
}


# Annotate CCGG/CG/SNP status
annotate_cg_status <- function(res, posfile, snptypefile){

    print("Loading pos table")
    # Add ccgg and cg status as two fields to the SNP loci
    pos = data.table::fread(posfile)

    # Format unique positions for annotation table
    print("Formatting SNP CG")
    SNPl = data.table()
    SNPl$chrom = res$Chromosome
    SNPl$start = res$Position
    SNPl$end = res$Position
    SNPl = unique(SNPl)

    # Annotate CG
    cgpos = pos[motif == 'CG', .(chrom, start, end)]; cgpos$is_cg=T; data.table::setkey(cgpos, chrom, start, end)
    ccggpos = pos[motif == 'CCGG',]; ccggpos$is_ccgg=T; data.table::setkey(ccggpos, chrom, start, end)

    # Annotate. For every item of SNPl, add first item of cgpos and ccggpos that overlaps
    print("Annotating SNP CG")
    cglo = data.table::foverlaps(SNPl, cgpos, type='any', which=TRUE, mult='first')
    SNPl$is_cg = ifelse(is.na(cglo), FALSE, TRUE)
    ccglo = data.table::foverlaps(SNPl, ccggpos, type='any', which=TRUE, mult='first')
    SNPl$is_ccgg = ifelse(is.na(ccglo), FALSE, TRUE)
    
    # Add to res 
    SNPl$Position = SNPl$start
    SNPl$Chromosome = SNPl$chrom
    res = merge(res, SNPl[, .(Chromosome, Position, is_cg, is_ccgg)], by=c("Chromosome", "Position"), all.x=T)

    # Load SNP types and add to res
    print("Annotating SNP type")
    snptype = data.table::fread(snptypefile)
    res = merge(res, snptype, by=c("Chromosome", "Position"), all.x=T)

    return(res)
}


# Write a function to get genotype score across samples
get_genotype_score <- function(f_wgs_normal_baf, f_bs_normal_baf) {
    stopifnot(length(f_wgs_normal_baf) == length(f_bs_normal_baf))
    
    results <- list()
    
    for (i in seq_along(f_wgs_normal_baf)) {


        # Load sample identifier
        sname = paste0('sample_', i)
        print(paste0("Processing scores for ", sname))

        # Load wgs and bsseq data
        wgs = data.table::fread(f_wgs_normal_baf[i])
        wgs$Chromosome = as.character(wgs$Chromosome)
        names(wgs)[3] = 'wgs_baf'
        bs = data.table::fread(f_bs_normal_baf[i])
        bs$Chromosome = as.character(bs$Chromosome)
        names(bs)[3] = 'bs_baf'
        
        # Merge datasets
        res = merge(wgs, bs, by = c('Chromosome', 'Position'))

        # Caclulate BAF difference
        res$bafdiff = res$wgs_baf - res$bs_baf
        
        # Add genotype metrics per SNP
        res$wgs_het = dplyr::between(res$wgs_baf, 0.1, 0.9)
        res$bs_het = dplyr::between(res$bs_baf, 0.1, 0.9)
        res[, bs_acc_01_09:= dplyr::case_when(
            wgs_het & bs_het ~ 'TP',
            !wgs_het & !bs_het ~ 'TN',
            wgs_het & !bs_het ~ 'FN',
            !wgs_het & bs_het ~ 'FP'
        )]
        
        # Add sample ID
        res$sample_id = sname
        
        results[[i]] <- res
    }
    
    return(data.table::rbindlist(results))
}

# Define allele_data for BAF from WGS allele frequency fieles
get_allele_data <- function(xv) {
  # Usage: allele_data = get_allele_data('data/meta/battenberg_hg38_allele_data.csv.gz')
  allele_data=data.table::fread(xv)
  names(allele_data)=c("#CHR", "POS", "ref_sel", "alt_sel")
  allele_data$POS= as.numeric(allele_data$POS)
  setkey(allele_data, "#CHR", "POS")
  return(allele_data)
}

# Get BAF from WGS allele frequency files
generate_unrandomised_battenberg_baf <- function(sample_id, allele_data) {
  normal_freq = rbindlist(lapply(Sys.glob(glue('/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_profiles/{sample_id}/A-GetAlleleCounts/S*_alleleFrequencies_chr*', sample_id=sample_id)), data.table::fread))
  setkey(normal_freq, "#CHR", "POS")
  
  norm = merge(normal_freq, allele_data, by=c("#CHR", "POS"), all.x=T)
  norm$ref_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$ref_sel)+2))])
  norm$alt_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$alt_sel)+2))])
  norm$BAF = norm$alt_count / (norm$ref_count + norm$alt_count)
  norm_table = norm[, .(Chromosome = `#CHR`, Position = POS, BAF)]
  norm_table = norm_table[!is.na(BAF)]
  names(norm_table)[3] = sample_id
  
  outdir = fs::dir_create(glue('out/wgs_baf/'))
  data.table::fwrite(norm_table, glue("{outdir}/{sample_id}_wgs_BAF.csv.gz"))
}

# Analysis of mCRPC SNP error ----

# Setup
# Specific to me (not lesky)
library(readr)

# Data
# Unfortunately, WGBS BAF has already been select for het SNPs and randomized.
# You can convert by taking the raw allele counts and filtering them:
## nbafdir = fs::dir_create('data/pre/240820_genotyping_accuracy')
## acdir='/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/archive/2023_data/230916_wgbs_ascat/wgbs_bb_tumor_normal/results/'
## bsraw <- fs::dir_ls(ac_dir, glob="*.N.SNPs.CpGs.all.sorted.csv.gz", recurse=T)
## bs_out = fs::path(nbafdir, paste0(str_extract(bsraw, 'DTB-\\d+.N'), '_normalBAF.tab'))
## for(i in seq_along(bsraw)){
##     allele_counts_to_battenberg_baf_tab(bsraw[i], bs_out[i])
## }


setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Genotyping_accuracy")
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
allele_data <- get_allele_data("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Genotyping_accuracy/battenerg_hg38_allele_data.csv.gz")
for (id in PT_ID){
    generate_unrandomised_battenberg_baf(id, allele_data)
}

nbafdir = fs::dir_create('out/wgbs_baf')
path_wgbs <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/06.CAMDAC/"
bsraw <- paste0(path_wgbs, PT_ID, "/Allelecounts/N/", PT_ID, ".N.SNPs.CpGs.all.sorted.csv.gz")
bs_out = paste0(nbafdir, "/", PT_ID, '_normalBAF.tab')
for(i in seq_along(bsraw)){
    allele_counts_to_battenberg_baf_tab(bsraw[i], bs_out[i])
}

# Normal BAF files
## WGS BAF files generated with scripts/scratch/240918_investigate_genotype_acc_err.R for mCRPC samples
f_wnbaf <- fs::dir_ls("out/wgs_baf/", glob = "*.csv.gz", recurse = T)
f_bnbaf <- fs::dir_ls("out/wgbs_baf/", glob = "*normalBAF.tab", recurse = T)

# Same order so same sample set
f_wnbaf <- f_wnbaf[order(f_wnbaf)]
f_bnbaf <- f_bnbaf[order(f_bnbaf)]

# Call the function to process scores (slow running)
genotype_scores <- get_genotype_score(f_wnbaf, f_bnbaf)
posfile='hg38_ccgg_cg_positions.csv.gz'
snptypefile='hg38_camdac_snp_ref_alt.csv.gz'
genotype_annotated <- annotate_cg_status(genotype_scores, posfile, snptypefile)
data.table::fwrite(genotype_annotated, 'genotype_scores_crc_wgbs.tsv', sep='\t')

# EDA of genotype accuracy ----
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(scales)
options(bitmapType='cairo')

# Load and format data
ga = data.table::fread('genotype_scores_crc_wgbs.tsv')
ga[, cg_status:=fcase(
    is_ccgg, 'CCGG',
    is_cg, 'CG',
    !is_ccgg & !is_cg, 'Non-CG'
)]
ga$bdiff_against_gs = ga$bafdiff # Cache default WGS - BSSEQ
ga$bafdiff = -1*ga$bafdiff # Flip as Lili is doing BSSEQ - WGS in her plots

# Generate plots
figdir = fs::dir_create('figures')
tpga = ga[bs_acc_01_09 %in% c("TP", "FN")]
tnga = ga[bs_acc_01_09 %in% c("TN", "FP")]

negs = ggplot(tnga, aes(x=bafdiff)) + geom_histogram() + facet_wrap(~cg_status, scales='free_y') + labs(x="WGBS BAF - WGS BAF", y="SNP counts (in m)") +
      scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 0.01, suffix = ""))
ggsave(path(figdir,'crc_genotype_accuracy_dist_tn.png'), negs, width=10, height=10)
poss  = ggplot(tpga, aes(x=bafdiff)) + geom_histogram() + facet_wrap(~cg_status, scales='free_y') + labs(x="WGBS BAF - WGS BAF", y="SNP counts (in m)") +
      scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 0.01, suffix = ""))
ggsave(path(figdir,'crc_genotype_accuracy_dist_tp.png'), poss, width=10, height=10)

# Calculate average FP and FN rates across all SNP types
fpr= ga[, .(
    FP_rate = sum(bs_acc_01_09 == 'FP') *100 / .N,
    FN_rate = sum(bs_acc_01_09 == 'FN') *100 / .N
), by=sample_id]
sprintf("Average FP rate: %f percent", mean(fpr$FP_rate))
sprintf("Average FN rate: %f percent", mean(fpr$FN_rate))
data.table::fwrite(fpr, path(figdir, 'crc_genotype_accuracy_fp_fn_rates.csv'))




# only plot for homozygous sites in WGS (BAF 0.15 ~ 0.85)
