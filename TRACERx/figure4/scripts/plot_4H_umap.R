#!/usr/bin/env Rscript
#
# Figure 4H: three UMAPs of one tumour's CpGs (epimutation class, m_t,1, m_t,2).
#
# PLOT ONLY. Loads precomputed per-CpG UMAP coordinates, clusters and classes from
# figure4/outputs/ASM_UMAP_recovered_from_pdf.RData (object `umap_recovered`); does
# not recompute UMAP / clustering / classification.

option_list <- list(
  optparse::make_option("--umap_rdata", type = "character",
    default = "~/data/CAMDAC_paper/figure4/outputs/ASM_UMAP_recovered_from_pdf.RData",
    help = "RData with object umap_recovered (per-patient CpG UMAP tables) [default %default]"),
  optparse::make_option("--tx_dir", type = "character", default = "~/rdeocd/txMeth"),
  optparse::make_option("--crukid", type = "character", default = "CRUK0084"),
  optparse::make_option("--work_dir", type = "character", default = "~/data/CAMDAC_paper/figure4")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })

load(path.expand(opt$umap_rdata))   # umap_recovered
if (!exists("umap_recovered")) stop("Expected object umap_recovered in ", opt$umap_rdata)
pi <- fread(file.path(path.expand(opt$tx_dir),
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))
patient <- unique(pi[CRUKID == opt$crukid, Patient])
stopifnot(length(patient) == 1, patient %in% names(umap_recovered))
df <- copy(umap_recovered[[patient]])
message(opt$crukid, " = ", patient, ": ", nrow(df), " CpGs")

class_lv <- c("clonal-AS-DMP1", "clonal-AS-DMP2", "clonal-bi-DMP",
              "subclonal-AS-DMP1", "subclonal-AS-DMP2", "subclonal-bi-DMP", "subclonal-ambiguous")
class_col <- c("clonal-AS-DMP1" = "purple3", "clonal-AS-DMP2" = "deeppink3",
               "clonal-bi-DMP" = "chartreuse4", "subclonal-AS-DMP1" = "mediumpurple2",
               "subclonal-AS-DMP2" = "lightpink2", "subclonal-bi-DMP" = "gold2",
               "subclonal-ambiguous" = "deepskyblue2")
class_lab <- c(expression("clonal AS-DMP"[1]), expression("clonal AS-DMP"[2]), "clonal bi-DMP",
               expression("subclonal AS-DMP"[1]), expression("subclonal AS-DMP"[2]),
               "subclonal bi-DMP", "unclassified")
df[, class := factor(as.character(class), levels = class_lv)]

lims <- list(x = range(df$umap1) + c(-0.4, 0.4), y = range(df$umap2) + c(-0.4, 0.4))
base <- theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
        axis.title = element_text(size = 11), plot.margin = margin(2, 4, 2, 2))
pt <- function(aes_col) geom_point(aes_col, alpha = 0.75, shape = 16, size = 1.1)
scx <- scale_x_continuous(limits = lims$x, expand = c(0, 0), name = "umap1")
scy <- scale_y_continuous(limits = lims$y, expand = c(0, 0), name = "umap2")

p_class <- ggplot(df, aes(umap1, umap2)) + pt(aes(col = class)) +
  scale_color_manual(name = "Epimutation class", values = class_col, labels = class_lab, drop = FALSE) +
  scx + scy + base + theme(legend.position = "none")
p_t1 <- ggplot(df, aes(umap1, umap2)) + pt(aes(col = m_t_1)) +
  scale_color_gradient(name = expression(m[t*","*1]), low = "lightblue", high = "purple3", limits = c(0, 1)) +
  scx + scy + base + theme(legend.position = "none", axis.title.y = element_blank())
p_t2 <- ggplot(df, aes(umap1, umap2)) + pt(aes(col = m_t_2)) +
  scale_color_gradient(name = expression(m[t*","*2]), low = "limegreen", high = "deeppink3", limits = c(0, 1)) +
  scx + scy + base + theme(legend.position = "none", axis.title.y = element_blank())

leg_theme <- theme(legend.position = "right", legend.text = element_text(size = 9.5),
                   legend.title = element_text(size = 10), legend.key.height = unit(0.9, "lines"),
                   legend.margin = margin(0, 0, 0, 0))
bar_guide <- guide_colourbar(direction = "horizontal", reverse = TRUE, title.position = "left",
                             barwidth = unit(6.2, "lines"), barheight = unit(0.8, "lines"))
leg_class <- get_legend(p_class + guides(color = guide_legend(override.aes = list(size = 2.6, alpha = 1))) + leg_theme)
leg_t1 <- get_legend(p_t1 + guides(color = bar_guide) +
                       theme(legend.position = "right", legend.text = element_text(size = 7),
                             legend.title = element_text(size = 10), legend.margin = margin(0, 0, 0, 0)))
leg_t2 <- get_legend(p_t2 + guides(color = bar_guide) +
                       theme(legend.position = "right", legend.text = element_text(size = 7),
                             legend.title = element_text(size = 10), legend.margin = margin(0, 0, 0, 0)))

panels <- plot_grid(p_class, p_t1, p_t2, nrow = 1, rel_widths = c(1.12, 1, 1), align = "h", axis = "tb")
legends <- plot_grid(leg_class, leg_t1, leg_t2, ncol = 1, rel_heights = c(2.4, 0.7, 0.7), align = "v")
title <- ggdraw() + draw_label(opt$crukid, x = 0.01, hjust = 0, size = 11, colour = "grey35")
fig <- plot_grid(title, plot_grid(panels, legends, nrow = 1, rel_widths = c(3.0, 1.3)),
                 ncol = 1, rel_heights = c(0.07, 1))

out <- file.path(path.expand(opt$work_dir), "outputs", paste0("Fig4H_UMAP_", opt$crukid, ".pdf"))
ggsave(out, fig, width = 8.6, height = 2.9)
ggsave(sub("\\.pdf$", ".png", out), fig, width = 8.6, height = 2.9, dpi = 200)
message("Wrote ", out)
