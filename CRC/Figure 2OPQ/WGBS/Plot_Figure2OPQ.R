# Figure 2OPQ updated: Phased WT reads compared with adjacent normal with LOH and on all copies using CAMDAC asm module ----------
# the input SNV file only contains LOH regions with clonal SNV and the mutation is filtered to be on all remaining copies.

# f_path <- "/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/"
# asm_paths <- Sys.glob(paste0(f_path, "07.phasing_asm/", PT_ID, "/Methylation/T/*.T.asm_meth.csv.gz"))
# t_paths <- Sys.glob(paste0(f_path, "06.CAMDAC/", PT_ID, "/Methylation/T/*.T.pure.csv.gz"))
# 
# phased_obj <- data.frame()
# for (i in 1:6){
#   dt_camdac <- data.table::fread(t_paths[i])
#   dt <- data.table::fread(asm_paths[i])
#   dt <- dt[, .(chrom, start, end, alt_m, ref_m, 
#                alt_cov = alt_total_counts_m, ref_cov = ref_total_counts_m)]
#   dt$chrom <- as.character(dt$chrom)
#   res <- inner_join(dt, dt_camdac, by = c("chrom", "start", "end"))
#   phased_obj <- rbind(phased_obj, res)
#   print(cor(res$ref_m[res$ref_cov >= 10], res$m_n[res$ref_cov >= 10], use = "complete.obs")) 
# }

# load saved object from above
phased_obj <- readRDS("Figure_2OPQ_WGBS_phased_obj.RDS")

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