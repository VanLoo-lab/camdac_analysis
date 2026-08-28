library(fs)
library(data.table)
library(stringr)
library(glue)
library(dplyr)
library(VariantAnnotation)

# Figure 2N: WGBS and WGS-derived VAF ----------

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
  tumor_freq = rbindlist(lapply(Sys.glob(glue('/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_profiles/{sample_id}/A-GetAlleleCounts/IPCT*_alleleFrequencies_chr*', sample_id=sample_id)), data.table::fread))
  setkey(tumor_freq, "#CHR", "POS")
  
  norm = merge(tumor_freq, allele_data, by=c("#CHR", "POS"), all.x=T)
  norm$ref_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$ref_sel)+2))])
  norm$alt_count = as.numeric(data.frame(norm)[cbind(1:nrow(norm), (as.numeric(norm$alt_sel)+2))])
  norm$BAF = norm$alt_count / (norm$ref_count + norm$alt_count)
  norm_table = norm[, .(Chromosome = `#CHR`, Position = POS, BAF)]
  norm_table = norm_table[!is.na(BAF)]
  names(norm_table)[3] = sample_id
  
  outdir = fs::dir_create(glue('wgs_vaf/'))
  data.table::fwrite(norm_table, glue("{outdir}/{sample_id}_wgs_VAF.csv.gz"))
}

# Convert allele counts to Battenberg *BAF table for WGBS
allele_counts_to_battenberg_baf_tab <- function(ac_file, outfile){
    # ac_file: A CAMDAC allele counts file
    # outfile: A file to write the Battenberg BAF table to
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

setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Fig2_NO")
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")

# Get WGS VAF for all the sites
allele_data <- get_allele_data("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Genotyping_accuracy/battenerg_hg38_allele_data.csv.gz")
for (id in PT_ID){
    generate_unrandomised_battenberg_baf(id, allele_data)
}

# Get WGBS VAF for all the sites
path_wgbs <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/06.CAMDAC/"
bsraw <- paste0(path_wgbs, PT_ID, "/Allelecounts/T/", PT_ID, ".T.SNPs.CpGs.all.sorted.csv.gz")
nbafdir = fs::dir_create('wgbs_baf')
bs_out = paste0(nbafdir, "/", PT_ID, '_tumorBAF.tab')
for(i in seq_along(bsraw)){
    allele_counts_to_battenberg_baf_tab(bsraw[i], bs_out[i])
}

# Only take SNVs
metadata <- read.csv("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/CAMDAC_vcf_metadata.csv")
PID <- metadata$PT_ID
VCF <- metadata$VCF2

vaf_compare <- data.table()
for(i in 1:6){
  vcf <- readVcf(VCF[i], genome = "hg38")
  is_snv <- nchar(ref(vcf)) == 1 & nchar(unlist(alt(vcf))) == 1
  vcf <- vcf[is_snv]
  tmp <- rowRanges(vcf)
  snv <- data.table(Chromosome = gsub("chr", "", as.character(seqnames(tmp))),
                    Position = start(tmp),
                    REF = as.character(mcols(tmp)$REF),
                    # ALT field has to be unlisted first
                    # there might be more than one alt allele
                    ALT = sapply(mcols(tmp)$ALT, function(x) as.character(x[1])))
  vaf_wgs <- data.table::fread(paste0("wgs_vaf/", PID[i], "_wgs_VAF.csv.gz"))
  vaf_snv <- merge(vaf_wgs, snv, all = FALSE)
  vaf_snv[, PT_ID := PID[i]]
  names(vaf_snv) <- c("Chromosome", "Position", "VAF_WGS", "REF", "ALT", "PT_ID")
  
  vaf_wgbs <- data.table::fread(paste0("wgbs_baf/", PID[i], "_tumorBAF.tab"))
  vaf_snv <- merge(vaf_wgbs, vaf_snv, all = FALSE)
  names(vaf_snv) <- c("Chr", "Position", "VAF_WGBS", "VAF_WGS", "REF", "ALT", "PT_ID")
  vaf_compare <- rbind(vaf_compare, vaf_snv)
}
data.table::fwrite(vaf_compare, "SNV_VAF_compare.csv.gz")
vaf_compare <- data.table::fread("SNV_VAF_compare.csv.gz")
cor(vaf_compare$VAF_WGS, vaf_compare$VAF_WGBS)

png("Figure_2N_VAF_comparison.png")
par(mar = c(5, 5, 3, 3))  # Reset outer margins and plot area
plot(vaf_compare$VAF_WGS, vaf_compare$VAF_WGBS,
    xlab = "WGS-derived VAF",
    ylab = "WGBS-derived VAF",
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 5606\nr = 0.86", side = 3, line = 0.5)
dev.off()


# add coverage 
dt_6 <- lapply(Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/06.CAMDAC/", PT_ID, "/Allelecounts/T/*.T.SNPs.CpGs.all.sorted.csv.gz")),
               function(file) {
                 tmp <- data.table::fread(file)
                 tmp <- tmp[, .(Chr = chrom, Position = start, total_counts, ref_counts, alt_counts)]
                 tmp[, PT_ID := strsplit(basename(file), "\\.")[[1]][1]]
               })
tmp <- do.call(rbind, dt_6)
vaf_compare$Chr <- as.character(vaf_compare$Chr)
tmp0 <- inner_join(tmp, vaf_compare)
cor(tmp0$VAF_WGS, tmp0$VAF_WGBS)

tmp_20 <- tmp0[tmp0$total_counts > 50, ]
cor(tmp_20$VAF_WGS, tmp_20$VAF_WGBS)
png("Figure_2N_VAF_comparison_cov50.png")
par(mar = c(5, 5, 3, 3))  # Reset outer margins and plot area
plot(tmp_20$VAF_WGS, tmp_20$VAF_WGBS,
    xlab = "WGS-derived VAF",
    ylab = "WGBS-derived VAF",
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 3585\nr = 0.89", side = 3, line = 0.5)
dev.off()


# Should take LOH regions only
vaf_compare <- data.table::fread("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Fig2_NO/SNV_VAF_compare.csv.gz")
loh_regions <- data.table()
for (i in 1:6) {
  cna <- data.table::fread(paste0(path_wgbs, "/", PT_ID[i], "/Copynumber/T/", PT_ID[i], ".T.cna.txt"))
  loh <- cna[cna$minor_cn == 0, 1:3]
  loh[, PT_ID := PT_ID[i]]
  loh_regions <- rbind(loh_regions, loh)
}
vaf_out <- vaf_compare[loh_regions,
  .(Chr, Position = x.Position,
    start = i.start,
    end = i.end, VAF_WGBS, VAF_WGS, PT_ID),
  on = .(Chr = chrom, PT_ID, Position >= start, Position <= end),
  nomatch = 0
]
data.table::fwrite(vaf_out, "SNV_VAF_compare_LOH.csv.gz")
vaf_out[, diff := VAF_WGBS - VAF_WGS]
cor(vaf_out$VAF_WGS, vaf_out$VAF_WGBS)
png("Figure_2N_VAF_comparison_LOH.png")
par(mar = c(5, 5, 3, 3))  # Reset outer margins and plot area
plot(vaf_out$VAF_WGS, vaf_out$VAF_WGBS,
    xlab = "WGS-derived VAF",
    ylab = "WGBS-derived VAF",
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 1123\nr = 0.92", side = 3, line = 0.5)
dev.off()

# Figure 2N RRBS --------------
# Only take SNVs
metadata <- read.csv("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/CAMDAC_vcf_metadata.csv")
PID <- metadata$PT_ID
VCF <- metadata$VCF2

vaf_compare <- data.table()
for(i in 1:6){
  vcf <- readVcf(VCF[i], genome = "hg38")
  is_snv <- nchar(ref(vcf)) == 1 & nchar(unlist(alt(vcf))) == 1
  vcf <- vcf[is_snv]
  tmp <- rowRanges(vcf)
  snv <- data.table(Chromosome = gsub("chr", "", as.character(seqnames(tmp))),
                    Position = start(tmp),
                    REF = as.character(mcols(tmp)$REF),
                    # ALT field has to be unlisted first
                    # there might be more than one alt allele
                    ALT = sapply(mcols(tmp)$ALT, function(x) as.character(x[1])))
  
  # WGS VAF
  vaf_wgs <- data.table::fread(paste0("wgs_vaf/", PID[i], "_wgs_VAF.csv.gz"))
  vaf_snv <- merge(vaf_wgs, snv, all = FALSE)
  vaf_snv[, PT_ID := PID[i]]
  names(vaf_snv) <- c("Chromosome", "Position", "VAF_WGS", "REF", "ALT", "PT_ID")
  
  # RRBS VAF
  load(Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/Batch*/bam_final/", 
    gsub("-0", "-", PID[i]), "/Copy_number/T1/*.T1.SNPs.RData")))
  vaf_rrbs <- dt[, c(1,2, 5, 6, 7)]
  names(vaf_rrbs) <- c("Chromosome", "Position", "REF", "ALT", "VAF_RRBS")

  # merge
  vaf_snv_rrbs <- merge(vaf_rrbs, vaf_snv, all = FALSE)
  names(vaf_snv_rrbs) <- c("Chr", "Position", "REF", "ALT", "VAF_RRBS", "VAF_WGS", "PT_ID")
  vaf_compare <- rbind(vaf_compare, vaf_snv_rrbs)
}
data.table::fwrite(vaf_compare, "SNV_VAF_compare_rrbs.csv.gz")

cor(vaf_compare$VAF_WGS, vaf_compare$VAF_RRBS)
png("Figure_2N_VAF_comparison_rrbs.png")
par(mar = c(5, 5, 3, 3))  # Reset outer margins and plot area
plot(vaf_compare$VAF_WGS, vaf_compare$VAF_RRBS,
    xlab = "WGS-derived VAF",
    ylab = "RRBS-derived VAF",
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 127\nr = 0.64", side = 3, line = 0.5)
dev.off()

vaf_filtered <- vaf_compare[vaf_compare$VAF_RRBS > 0.1, ]
cor(vaf_filtered$VAF_WGS, vaf_filtered$VAF_RRBS)
png("Figure_2N_VAF_comparison_rrbs_minvaf0.1.png")
par(mar = c(5, 5, 3, 3))  # Reset outer margins and plot area
plot(vaf_filtered$VAF_WGS, vaf_filtered$VAF_RRBS,
    xlab = "WGS-derived VAF",
    ylab = "RRBS-derived VAF",
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 47\nr = 0.83", side = 3, line = 0.5)
dev.off()



# check why RRBS has so few SNVs with VAF
# LORC-02
# SNV: 18932
# overlap of vaf_wgs and snv: 1247
# WGBS SNP sites: 27120821
# WGBS VAF: 1073
# RRBS SNP sites: 2353448
# RRBS VAF: 31
# The ratio here is a little strange: should have more SNVs in CGIs?
# Lung epitracerx data: Figure 2N showed ~ 36 SNVs per sample (are these phased SNVs?) 
# (SNVs are filtered with coverage > 30 and ALT > 5; maybe also LOH region only as well?)
# Original data ~ 120 phased SNVs per lung sample; CRC has 25% mutation rate so ~40 phased SNV per patient without filtering

# Figure 2O: Phased WT reads compared with adjacent normal ----------
f_path <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/"
wt_paths <- Sys.glob(paste0(f_path, "07.phasing/", PT_ID, "/*/Methylation/WT/*.WT.m.csv.gz"))
snv_paths <- Sys.glob(paste0(f_path, "07.phasing/", PT_ID, "/*/Methylation/SNV/*.SNV.m.csv.gz"))
n_paths <- Sys.glob(paste0(f_path, "06.CAMDAC/", PT_ID, "/Methylation/N/*.N.m.csv.gz"))
t_paths <- Sys.glob(paste0(f_path, "06.CAMDAC/", PT_ID, "/Methylation/T/*.T.pure.csv.gz"))

phased_wt_obj <- phased_snv_obj <- data.frame()
for (i in 1:6){
  m_snv <- data.table::fread(snv_paths[i])
  m_snv$chrom <- as.character(m_snv$chrom)
  m_t <- data.table::fread(t_paths[i])
  phased_snv <- inner_join(m_snv, m_t, by = c("chrom", "start", "end"))
  phased_snv <- phased_snv[phased_snv$cov.x > 10 & phased_snv$nB == 0, ]
  phased_snv_obj <- rbind(phased_snv_obj, phased_snv)

  m_wt <- data.table::fread(wt_paths[i])
  m_wt$chrom <- as.character(m_wt$chrom)
  m_n <- data.table::fread(n_paths[i])
  phased_wt <- inner_join(m_wt, m_n, by = c("chrom", "start", "end"))
  phased_wt <- inner_join(phased_wt, m_t[, c("chrom", "start", "end", "nB")])
  phased_wt <- phased_wt[phased_wt$cov.x > 10 & phased_wt$nB == 0, ]
  phased_wt_obj <- rbind(phased_wt_obj, phased_wt)
}


plot(phased_wt_obj$m.x, phased_wt_obj$m.y,
    xlab =expression(paste("WT m"["n"],"")),
    ylab = expression(paste("m"["n"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))



plot(phased_snv_obj$m.x, phased_snv_obj$m_t,
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["t"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))

plot(phased_snv_obj$m.x, phased_snv_obj$m.y,
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["b"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))


# Figure 2OPQ updated: Phased WT reads compared with adjacent normal with LOH and on all copies using CAMDAC asm module ----------
# this should be the final version for WGBS figure 2OPQ
setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Fig2_NO")
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
f_path <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/"
asm_paths <- Sys.glob(paste0(f_path, "07.phasing_asm/", PT_ID, "/Methylation/T/*.T.asm_meth.csv.gz"))
t_paths <- Sys.glob(paste0(f_path, "06.CAMDAC/", PT_ID, "/Methylation/T/*.T.pure.csv.gz"))

phased_obj <- data.frame()
for (i in 1:6){
  dt_camdac <- data.table::fread(t_paths[i])
  dt <- data.table::fread(asm_paths[i])
  dt <- dt[, .(chrom, start, end, alt_m, ref_m, 
               alt_cov = alt_total_counts_m, ref_cov = ref_total_counts_m)]
  dt$chrom <- as.character(dt$chrom)
  res <- inner_join(dt, dt_camdac, by = c("chrom", "start", "end"))
  phased_obj <- rbind(phased_obj, res)
  print(cor(res$ref_m[res$ref_cov >= 10], res$m_n[res$ref_cov >= 10], use = "complete.obs")) 
}
saveRDS(phased_obj, file = "Figure_2OPQ_WGBS_phased_obj.RDS")

# m_wt vs m_n
pdf("Figure_2O_WGBS_asmmodule.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(phased_obj$ref_m[phased_obj$ref_cov >= 20], phased_obj$m_n[phased_obj$ref_cov >= 20],
    xlab =expression(paste("WT m"["n"],"")),
    ylab = expression(paste("m"["n"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
abline(0,1, lwd=5, lty=2, col="red")
mtext("n = 4607\nr = 0.82", side = 3, line = 0.5)
dev.off()
cor(phased_obj$ref_m[phased_obj$ref_cov >= 10], phased_obj$m_n[phased_obj$ref_cov >= 10], use = "complete.obs")

# m_snv vs m_decon
pdf("Figure_2P_WGBS_asmmodule.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(phased_obj$alt_m[phased_obj$alt_cov >= 10], phased_obj$m_t[phased_obj$alt_cov >= 10],
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["t"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
abline(0,1, lwd=5, lty=2, col="red")
dev.off()
cor(phased_obj$alt_m[phased_obj$alt_cov >= 10], phased_obj$m_t[phased_obj$alt_cov >= 10], use = "complete.obs")

# m_snv vs m_bulk
pdf("Figure_2Q_WGBS_asmmodule.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(phased_obj$alt_m[phased_obj$alt_cov >= 10], phased_obj$m[phased_obj$alt_cov >= 10],
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["b"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
abline(0,1, lwd=5, lty=2, col="red")
dev.off()
cor(phased_obj$alt_m[phased_obj$alt_cov >= 10], phased_obj$m[phased_obj$alt_cov >= 10], use = "complete.obs")



# WGBS Figure 2N with phasing objects --------
# this is the correct version for WGBS Fig 2N
setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Fig2_NO")
library(tidyverse)
library(data.table)
metainfo <- read_csv("/rsrch6/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/CAMDAC_vcf_metadata.csv")
res <- data.table()
for (i in 1:6){
  pt <- metainfo$PT_ID[i]
  fpath <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/07.phasing_asm_snv/"
  df1 <- fread(paste0(fpath, pt, "/AlleleSpecific/T/", pt, ".T.asm_hap_stats.csv.gz"))
  df1 <- df1[, 1:7]
  # hap_POS and hap_id corresponds to T.asm_snps.txt in the same folder, which are actually SNVs now but formated as SNPs to run phasing
  colnames(df1) <- c("chrom", "ref", "alt", "pos", "hap_id", "BAF_wgbs", "hap_reads")
  df1$chrom <- paste0("chr", df1$chrom)
  df2 <- fread(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/vcf_snv/", pt, "_vaf1.txt"))
  df <- inner_join(df1, df2)
  print(nrow(df))
  df$PT_ID <- pt
  res <- rbind(res, df)
}

res0 <- res[res$hap_reads > 30 &
            res$alt_count > 5 &
            res$hap_reads * res$BAF_wgbs > 5 &
            res$ref_count + res$alt_count > 30]
saveRDS(res, file = "Figure_2N_WGBS.RDS")
saveRDS(res0, file = "Figure_2N_WGBS_cov30_alt5.RDS")

pdf("Figure_2N_phasing.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res$BAF2, res$BAF_wgbs,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("WGBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 88832\nr = 0.84", side = 3, line = 0.5)
dev.off()

pdf("Figure_2N_phasing_cov30_alt5.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res0$BAF2, res0$BAF_wgbs,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("WGBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 47887\nr = 0.83", side = 3, line = 0.5)
dev.off()

cor(res$BAF2, res$BAF_wgbs)
cor(res0$BAF2, res0$BAF_wgbs)

# Filter by LOH, clonal and MCN
res0 <- readRDS("Figure_2N_WGBS_cov30_alt5.RDS")
tmp <- rbindlist(lapply(Sys.glob("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/vcf_clonal/*_LOH_clonal_allcopy_vaf.txt"), function (file) {
  pt <- sub("_LOH_clonal_allcopy_vaf.txt$", "", basename(file))
  data <- fread(file)
  data[, PT_ID := pt]  # Add the sample ID as a new column
  return(data)
}))
res0 <- inner_join(res0, tmp)
saveRDS(res0, file = "Figure_2N_WGBS_clonal_LOH_allcopy_cov30_alt5.RDS")

pdf("Figure_2N_phasing_clonal_LOH_allcopy_cov30_alt5.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res0$BAF2, res0$BAF_wgbs,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("WGBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 2983\nr = 0.91", side = 3, line = 0.5)
dev.off()

cor(res0$BAF2, res0$BAF_wgbs)


# check number of phased CpGs (# phased SNVs is the count values in the plots above)
nsnv <- ncpg <- 0
res <- NULL
for (i in 1:6){
  pt <- metainfo$PT_ID[i]
  fpath <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/"
  dt1 <- fread(paste0(fpath, "07.phasing_asm_snv/", pt, "/Methylation/T/", pt, ".T.asm_meth.csv.gz"))
  dt2 <- fread(paste0(fpath, "07.phasing_asm_snv/", pt, "/AlleleSpecific/T/", pt, ".T.asm_phase_map.csv.gz"))
  dt2 <- unique(dt2[, 2:5])
  dt3 <- fread(paste0(fpath, "07.phasing_asm_snv/", pt, "/AlleleSpecific/T/", pt, ".T.asm_hap_stats.csv.gz"))
  dt <- inner_join(dt1, dt2) # could have CpGs that are phased to multiple SNVs
  dt$width <- dt$end - dt$start + 1
  dt <- inner_join(dt, dt3)
  dt$PT_ID <- pt
  
  # add info from bulk deconvolution
  dt_bulk <- fread(paste0(fpath, "06.CAMDAC/", pt, "/Methylation/T/", pt, ".T.pure.csv.gz"))
  purity <- dt_bulk$purity[1]
  dt_tmp <- dt_bulk[, .(chrom, start, end, nA, nB, M_n, UM_n, m_n, m_b = m, m_t)]                 
  dt_all <- as.data.table(inner_join(dt, dt_tmp))
  dt_phased <- dt_all[order(chrom, start, end, width, nA, nB),
            .(PT_ID, genotype_1 = "ref", genotype_2 = "alt", POS_genotype = hap_POS, 
              M_1 = ref_total_counts_m * ref_m, UM_1 = ref_total_counts_m * (1 - ref_m), total_counts_m_1 = ref_total_counts_m, 
              M_2 = ref_total_counts_m * ref_m, UM_2 = ref_total_counts_m * (1 - ref_m), total_counts_m_2 = alt_total_counts_m,
              M_n, UM_n, 
              m_1 = ref_m, m_2 = alt_m, m_n, m_b, m_t),
              keyby=c("chrom","start","end","width","nA", "nB")]
  dt_phased[, SNV_VAF_methylation := as.numeric(total_counts_m_2)/(as.numeric(total_counts_m_2)+as.numeric(total_counts_m_1))]
  dt_phased[, SNV_mcn_methylation := (1/purity)*SNV_VAF_methylation*((nA+nB)*purity+2*(1-purity))] 
  dt_phased <- dt_phased[!is.na(SNV_VAF_methylation) & !is.na(SNV_mcn_methylation), ]

  # add info from WGS
  wgs <- fread(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/vcf_snv/", pt, "_vaf1.txt"))
  wgs1 <- fread(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/vcf_clonal/", pt, "_clonal_vaf.txt"))
  wgs <- inner_join(wgs, wgs1[, 1:4])
  
  tmp <- wgs[, .(chrom = gsub("chr", "", chrom), POS_genotype = pos, SNV_VAF_WGS = BAF2, ref_count, alt_count)]
  # tmp <- tmp[alt_count > 3 & ref_count + alt_count > 18 & SNV_VAF_WGS > 0.1,] # WGS has lower coverage comapre to lili's data, adjust threshold
  tmp <- tmp[alt_count > 3 & SNV_VAF_WGS > 0.1, ] # WGS has lower coverage comapre to lili's data, adjust threshold
  dt_final <- as.data.table(inner_join(dt_phased, tmp[, 1:3]))
  dt_final[, SNV_mcn := (1/purity)*SNV_VAF_WGS*((nA+nB)*purity+2*(1-purity))]
  dt_final_clonal_loh_mcn <- dt_final[nB==0 & nA-0.3 <= SNV_mcn & nA+0.3 > SNV_mcn, ]
  res <- rbind(res, dt_final_clonal_loh_mcn)
  ## res <- rbind(res, dt_final)
}

## saveRDS(res, file = "Figure_2N_WGBS_no_filtering.RDS")

saveRDS(res, file = "Figure_2N_WGBS_final_VAF0.1_ALT3.RDS")
nrow(unique(res[, .(PT_ID, POS_genotype)])) # number of phased SNVs
nrow(res) # number of phased CpGs
cor(res$SNV_VAF_WGS, res$SNV_VAF_methylation)

# pdf("Figure_2N_phasing_no_filter.pdf")
# par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
# plot(res$SNV_VAF_WGS, res$SNV_VAF_methylation,  xlab =expression(paste("WGS-derived VAF")),
#     ylab = expression(paste("WGBS-derived VAF")),
#     xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
#     las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
#     col=rgb(0, 0, 0, alpha = 0.5))
# mtext("n = 65894\nr = 0.7632202", side = 3, line = 0.5)
# dev.off()

pdf("Figure_2N_phasing_clonal_LOH_allcopy_final.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res$SNV_VAF_WGS, res$SNV_VAF_methylation,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("WGBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 2779\nr = 0.88", side = 3, line = 0.5)
dev.off()

res0 <- res[total_counts_m_2 > 5 & total_counts_m_2 + total_counts_m_1 > 30, ]
nrow(unique(res0[, .(PT_ID, POS_genotype)])) # number of phased SNVs
nrow(res0) # number of phased CpGs
cor(res0$SNV_VAF_WGS, res0$SNV_VAF_methylation)

pdf("Figure_2N_phasing_clonal_LOH_allcopy_final_cov30_alt5.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res0$SNV_VAF_WGS, res0$SNV_VAF_methylation,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("WGBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 1648\nr = 0.92", side = 3, line = 0.5)
dev.off()

# res00 <- res0[res$hap_reads > 30 &
#             res$alt_count > 5 &
#             res$hap_reads * res$BAF_wgbs > 5 &
#             res$ref_count + res$alt_count > 30]
# pdf("Figure_2N_phasing_clonal_LOH_allcopy_cov_30_alt_5_redo.pdf")
# par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
# plot(res00$SNV_VAF_WGS, res00$hap_BAF,  xlab =expression(paste("WGS-derived VAF")),
#     ylab = expression(paste("WGBS-derived VAF")),
#     xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
#     las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 0.5,
#     col=rgb(0, 0, 0, alpha = 0.5))
# mtext("n = 2733\nr = 0.89", side = 3, line = 0.5)
# dev.off()
# cols <- c("PT_ID", "chrom", "start", "POS_genotype", "width",
#         "m_n", "M_n", "UM_n",
#         "m_b", "m_t",
#         "m_1", "m_2", "M_2", "UM_2",
#         "m_t_2",
#         "SNV_mcn_methylation",
#         "SNV_VAF_WGS", "SNV_VAF_methylation")
# cols_old <- c("PT_ID", "chrom", "start", "hap_POS", "width",
#         "m_n", "M_n", "UM_n",
#         "m_b", "m_t",
#         "alt_m", "ref_m", "M_2", "UM_2",
#         "m_t_2",
#         "SNV_mcn_methylation",
#         "SNV_VAF_WGS", "SNV_VAF_methylation")


# Figure 2NO for RRBS using phased SNVs -----------
# load snv_phased rrbs object (already filtered by clonal, LOH and MCN with Lili's phasing code)
df <- readRDS("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/RRBS_phasing/phased_objs/df_phased_LORC-2.rds")
df_tumour <- df[!is.na(total_counts_m), ]
df_phased_methylation <- df_tumour[order(chrom,start,end,width,POS,genotype),
                                        .(genotype_1="ref",genotype_2="alt", num_genotypes=length(genotype),
                                        POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        SNV=SNP[genotype=="ref"&SNP%in%SNP[genotype=="alt"]],
                                        counts_1=total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        counts_2=total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                        M_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                        UM_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        total_counts_m_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        M_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                        UM_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                        total_counts_m_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]),
                                        keyby=c("CHR","chrom","start","end","width","POS","ref","alt")]

df_phased_methylation[, SNV_VAF_methylation := as.numeric(total_counts_m_2)/
                        (as.numeric(total_counts_m_2)+as.numeric(total_counts_m_1))]

ids <- c("LORC-2", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
pts <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
path_phased <- "/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/RRBS_phasing/"
df_vaf <- data.frame()
for (i in 1:6){
  df <- readRDS(paste0(path_phased, "df_phased_", ids[i], ".rds"))
  df_tumour <- df[!is.na(total_counts_m), ]
  df_phased_methylation <- df_tumour[order(chrom,start,end,width,POS,genotype),
                                        .(genotype_1="ref",genotype_2="alt", num_genotypes=length(genotype),
                                        POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        SNV=SNP[genotype=="ref"&SNP%in%SNP[genotype=="alt"]],
                                        counts_1=total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        counts_2=total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                        M_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                        UM_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        total_counts_m_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                        M_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                        UM_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                        total_counts_m_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]),
                                        keyby=c("CHR","chrom","start","end","width","POS","ref","alt")]
  df_phased_methylation[, m_1 := as.numeric(M_1)/as.numeric(total_counts_m_1)] 
  df_phased_methylation[, m_2 := as.numeric(M_2)/as.numeric(total_counts_m_2)]
  df_phased_methylation[, SNV_VAF_methylation := as.numeric(total_counts_m_2)/
                        (as.numeric(total_counts_m_2)+as.numeric(total_counts_m_1))]
  df_phased_methylation$PT_ID <- pts[i]             
  df_vaf <- rbind(df_vaf, df_phased_methylation)
}

# Fig 2N
df_vaf1 <- df_vaf[, .(chrom, pos=as.integer(POS_genotype), total_counts_m_1, total_counts_m_2, SNV_VAF_methylation, PT_ID)]
# many of the dots here are not unique since more than 1 CpG could be phased to a SNV
df_vaf1$chrom <- paste0("chr", df_vaf1$chrom)
# add WGS VAF
tmp <- rbindlist(lapply(Sys.glob("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/vcf_clonal/*_LOH_clonal_allcopy_vaf.txt"), function (file) {
  pt <- sub("_LOH_clonal_allcopy_vaf.txt$", "", basename(file))
  data <- fread(file)
  data[, PT_ID := pt]  # Add the sample ID as a new column
  return(data)
}))

# my filtering rules of clonal, LOH and MCN further removed some sites
df_vaf2 <- inner_join(df_vaf1, tmp)
setnames(df_vaf2, old = "BAF", new = "VAF_WGS")
cor(df_vaf2$VAF_WGS, df_vaf2$SNV_VAF_methylation)

pdf("Figure_2N_phasing_clonal_LOH_allcopy_RRBS.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(df_vaf2$VAF_WGS, df_vaf2$SNV_VAF_methylation,  xlab =expression(paste("WGS-derived VAF")),
    ylab = expression(paste("RRBS-derived VAF")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1,
    col=rgb(0, 0, 0, alpha = 0.5))
mtext("n = 149\nr = 0.96", side = 3, line = 0.5)
dev.off()

saveRDS(df_vaf2, file = "Figure_2N_RRBS.RDS")

# Fig 2O
# add matched normal
df_vaf3 <- df_vaf[, .(chrom = CHR, start, end, genotype_1, genotype_2, m_1, m_2, total_counts_m_1, total_counts_m_2, PT_ID)]
res <- data.table()
for (i in 1:6){
  load(Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/Batch*/bam_final/", ids[i], "/Methylation/N1/dt_normal_m.RData")))
  tmp <- dt_normal_m
  tmp$PT_ID <- pts[i]
  tmp <- tmp[, .(chrom = CHR, start, end, m_n, cov_n)]
  df <- inner_join(tmp, df_vaf3)
  res <- rbindlist(list(res, df))
}

res <- res[!is.na(m_1), ]
res <- res[!is.na(m_2), ]
cor(res$m_1, res$m_n)
pdf("Figure_2O_RRBS.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res$m_1, res$m_n,
    xlab =expression(paste("WT m"["n"],"")),
    ylab = expression(paste("m"["n"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
abline(0,1, lwd=5, lty=2, col="red")
mtext("n = 1297\nr = 0.97", side = 3, line = 0.5)
dev.off()
saveRDS(res, file = "Figure_2O_RRBS.RDS")

res <- res[res$total_counts_m_1 > 30 & res$cov_n > 30, ]
cor(res$m_1, res$m_n)
pdf("Figure_2O_RRBS_cov30.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(res$m_1, res$m_n,
    xlab =expression(paste("WT m"["n"],"")),
    ylab = expression(paste("m"["n"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col=rgb(0, 0, 0, alpha = 0.5))
abline(0,1, lwd=5, lty=2, col="red")
mtext("n = 896\nr = 0.98", side = 3, line = 0.5)
dev.off()

saveRDS(res, file = "Figure_2O_RRBS_cov30.RDS")
View(res)


# Figure 2PQ for RRBS with updated phasing objects with all steps ------------
f_paths <- "/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/RRBS_phasing/phased_final/"
ids <- c("LORC-2", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
pts <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
dt <- data.table()
for (i in 1:6){
  tmp <- readRDS(paste0(f_paths, "asm_", ids[i], ".rds"))
  tmp$PT_ID <- pts[i]
  dt <- rbind(dt, tmp)
}
dt <- dt[dt$m_diff_tn2 > 0.2, ]


group_colors <- c("LORC-02" = "red", "LORC-60" = "orange", "LORC-62" = "yellow",
                  "LORC-68" = "green", "YORC-25" = "blue", "YORC-41" = "purple")
pdf("Figure_2P_RRBS_final.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(dt$m_t_2, dt$m_t,
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["t"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col = group_colors[dt$PT_ID])
abline(0,1, lwd=5, lty=2, col="red")
dev.off()

pdf("Figure_2Q_RRBS_final.pdf")
par(mar = c(4.1, 5, 4.1, 1.9), xaxs = "i", yaxs = "i", lwd=5)
plot(dt$m_t_2, dt$m_b,
    xlab =expression(paste("SNV m"["t"],"")),
    ylab = expression(paste("CAMDAC m"["b"],"")),
    xlim = c(0, 1.02), ylim = c(0, 1.02), bty = "l",
    las = 1, cex.axis = 1, pch=16, cex.lab=2, cex.axis=1.5, cex = 1.5,
    col = group_colors[dt$PT_ID])
abline(0,1, lwd=5, lty=2, col="red")
dev.off()
saveRDS(dt, file = "Figure_2PQ_RRBS.rds")