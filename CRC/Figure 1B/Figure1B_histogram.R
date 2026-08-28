# Histogram of bulk methylation rate by purity and CN states

library(tidyverse)
library(fst)
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
wgbs_paths <- Sys.glob(paste0("/rsrch6/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/T/", PT_ID, ".T.CAMDAC_results_per_CpG.fst"))

# Purity = 0.81, ploidy = 3.23
dt <- read.fst(wgbs_paths[3], as.data.table = TRUE)
# CpGs confidently unmethylated in the adjacent normal sample
dt <- dt[dt$m_x_high_n <= 0.2 & dt$m_x_low_n > 0, ]
CN_2 <- dt[dt$nA == 2 & dt$nB == 0, ]
CN_3 <- dt[dt$nA == 2 & dt$nB == 1, ]
CN_4 <- dt[dt$nA == 2 & dt$nB == 2, ]


png("Fig1B/LORC-62.png", width = 800, height = 1600)
par(mfrow = c(3, 1))
hist(CN_2$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 600))
hist(CN_3$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 5000))
hist(CN_4$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 300))
dev.off()
par(mfrow = c(1, 1))
save(CN_2, CN_3, CN_4, file = "Fig1B/LORC-62_3CN.RData")


# Purity = 0.64, ploidy = 4.52
dt <- read.fst(wgbs_paths[1], as.data.table = TRUE)
# CpGs confidently unmethylated in the adjacent normal sample
dt <- dt[dt$m_x_high_n <= 0.2 & dt$m_x_low_n > 0, ]
CN_2 <- dt[dt$nA == 2 & dt$nB == 0, ]
CN_3 <- dt[dt$nA == 2 & dt$nB == 1, ]
CN_4 <- dt[dt$nA == 2 & dt$nB == 2, ]


png("Fig1B/LORC-02.png", width = 800, height = 1600)
par(mfrow = c(3, 1))
hist(CN_2$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 200))
hist(CN_3$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 500))
hist(CN_4$m, breaks = 40, xlim = c(0, 1), ylim = c(0, 10000))
dev.off()
par(mfrow = c(1, 1))
save(CN_2, CN_3, CN_4, file = "Fig1B/LORC-02_3CN.RData")
