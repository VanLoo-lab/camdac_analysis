
# WGBS Figure 2N with phasing objects --------
setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Fig2_NO")
metainfo <- read_csv("/rsrch6/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/CAMDAC_vcf_metadata.csv")
res <- data.table()
for (i in 1:6){
  pt <- metainfo$PT_ID[i]
  fpath <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/07.phasing_asm_snv/"
  df1 <- fread(paste0(fpath, pt, "/AlleleSpecific/T/", pt, ".T.asm_hap_stats.csv.gz"))
  df1 <- df1[, 1:7]
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