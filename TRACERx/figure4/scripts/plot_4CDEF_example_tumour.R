#!/usr/bin/env Rscript
#
# Figure 4C-F: allele-specific methylation of one example tumour (CRUK0084 = CRUK0084).
#
#   4C  histogram of CpG counts vs m_t (R1, ASCN 1_1, hyper-methylated CpGs)
#   4D  m_t,1 vs m_t,2 (R1, same CpGs)
#   4E  m_t of R1 vs R3 / R1 vs R2 (CpGs shared across regions = the UMAP CpG set of Fig 4H)
#   4F  m_t,1 - m_t,2 of R1 vs R3 / R1 vs R2 (same CpGs)
#
# PLOT ONLY: reads the per-region DMP_classification_<sample>.RData written by
# compute_DMP_classification.R (unchanged) and the recovered Jul-2021 UMAP CpG set
# (outputs/ASM_UMAP_recovered_from_pdf.RData) to select the 4E/4F CpGs. Nothing is overwritten.
#
# Source: Notes_ASM_plots_14Jul2021.R (4C/D: L160-290 sample plots; 4E/F: L755-905).
# Circles / arrows / text are hand annotations in the paper; they are re-drawn here in data
# coordinates (switch off with --annotate FALSE). The chromosome cartoons in 4C are not drawn.
#
# NB: the Jul-2021 Notes code assigned the 4F axes from swapped columns (`m_t_1 := m_t_2_<sample>`),
# so the published panel actually showed m_t,2 - m_t,1 under a "m_t,1 - m_t,2" label. Here the axes are
# computed as labelled (m_t,1 - m_t,2 = m_t_1_raw - m_t_2_raw), so the AS-DMP1 clusters (m_t,1 > m_t,2)
# are top-right and AS-DMP2 bottom-left, i.e. the colours of the circles are the mirror image of the
# 2021 panel.

## Prefer CRUKID in defaults; resolve to Patient (LTX*) via clinical table for paths
.shared <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", "_shared", "resolve_tx_ids.R")
.shared <- normalizePath(.shared, mustWork = FALSE)
if (!file.exists(.shared)) {
  .shared <- normalizePath(file.path("..", "..", "_shared", "resolve_tx_ids.R"), mustWork = FALSE)
}
if (!exists("resolve_tx_patient_id", mode = "function")) {
  if (file.exists(.shared)) source(.shared) else stop("Missing TRACERx/_shared/resolve_tx_ids.R")
}
.tx_map <- load_tx_patient_map()
option_list <- list(
  optparse::make_option("--tx_dir", type = "character", default = "~/rdeocd/txMeth"),
  optparse::make_option("--work_dir", type = "character", default = "~/data/CAMDAC_paper/figure4"),
  optparse::make_option("--recovered_rdata", type = "character",
    default = "~/data/CAMDAC_paper/figure4/outputs/ASM_UMAP_recovered_from_pdf.RData"),
  optparse::make_option("--crukid", type = "character", default = "CRUK0084"),
  optparse::make_option("--ascn", type = "character", default = "1_1"),
  optparse::make_option("--annotate", type = "logical", default = TRUE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })

tx_dir <- path.expand(opt$tx_dir)
out_dir <- file.path(path.expand(opt$work_dir), "outputs")
clin <- fread(file.path(tx_dir, "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))
patient <- unique(clin[CRUKID == opt$crukid, Patient]); stopifnot(length(patient) == 1)
samples <- clin[Patient == patient, Sample][order(clin[Patient == patient, cov_t], decreasing = TRUE)]
message(opt$crukid, " = ", patient, "; samples: ", paste(samples, collapse = ", "))

## ---- per-region data (Notes L62-120) -------------------------------------------------
read_region <- function(s) {
  e <- new.env()
  load(file.path(tx_dir, patient, "Allele_specific_analysis", s, paste0("DMP_classification_", s, ".RData")), envir = e)
  dt <- as.data.table(e$dt_DMP_classification)
  dt <- dt[n1 != 0 & n2 != 0 & m_n_high < 0.1]                       # hyper-methylated in normal
  dt[, `:=`(chrom = as.character(chrom),
            m_t = (m_t_1_raw * n1 + m_t_2_raw * n2) / (n1 + n2),
            m_diff = m_t_1_raw - m_t_2_raw)]
  dt[, .(chrom, start, end, multip, m_t, m_t_1 = m_t_1_raw, m_t_2 = m_t_2_raw, m_diff)]
}
reg <- setNames(lapply(samples, read_region), samples)

s1 <- samples[1]
cd <- reg[[s1]][multip == opt$ascn]
message("4C/D: ", s1, " multip ", opt$ascn, ": ", nrow(cd), " CpGs")

## ---- shared CpG set for 4E/F ----------------------------------------------------------
load(path.expand(opt$recovered_rdata))                                # umap_recovered
keys <- unique(umap_recovered[[patient]][, .(chrom = as.character(chrom), start, end)])
ef <- Reduce(function(a, b) merge(a, b, by = c("chrom", "start", "end")),
             lapply(samples, function(s) {
               x <- merge(keys, reg[[s]], by = c("chrom", "start", "end"))
               x[, .(chrom, start, end, m_t, m_diff)][, setnames(.SD, c("m_t", "m_diff"), paste0(c("m_t_", "m_diff_"), s))]
             }))
message("4E/F: ", nrow(ef), " shared CpGs (UMAP set: ", nrow(keys), ")")
for (s in samples) set(ef, j = paste0("f_", s), value = ef[[paste0("m_diff_", s)]])   # m_t,1 - m_t,2

## ---- styling ---------------------------------------------------------------------------
col_pt <- "grey25"
cols <- c(mag = "deeppink3", pur = "purple3", grn = "chartreuse4", pink = "lightpink2", lav = "mediumpurple2")
ell <- function(cx, cy, a, b = a, col, n = 120, lw = 0.9) {
  th <- seq(0, 2 * pi, length.out = n)
  geom_path(data = data.frame(x = cx + a * cos(th), y = cy + b * sin(th)), aes(x, y), colour = col,
            linewidth = lw, inherit.aes = FALSE)
}
arr <- function(x, y, xend, yend) {
  geom_segment(data = data.frame(x, y, xend, yend), aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.09, "in"), type = "closed"), linewidth = 0.4,
               colour = "black", inherit.aes = FALSE)
}
txt <- function(x, y, label, col = "black", size = 3.2, hjust = 0.5, fontface = "plain", lineheight = 0.95) {
  annotate("text", x = x, y = y, label = label, colour = col, size = size, hjust = hjust,
           fontface = fontface, lineheight = lineheight)
}
ann <- function(...) if (opt$annotate) list(...) else NULL

theme_c <- theme_classic(base_size = 11) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 9), plot.margin = margin(14, 8, 4, 4))
theme_g <- theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
        panel.grid.minor = element_line(colour = "grey93", linewidth = 0.3),
        axis.line = element_line(colour = "black", linewidth = 0.5), axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 9), aspect.ratio = 1, plot.margin = margin(4, 6, 4, 4))
brk <- seq(0, 1, 0.2)
pt <- function(...) geom_point(alpha = 0.5, size = 1.1, shape = 16, colour = col_pt, ...)

## ---- 4C ---------------------------------------------------------------------------------
ymax_c <- 115
p4C <- ggplot(cd, aes(m_t)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "grey75", colour = "grey30", linewidth = 0.3) +
  scale_x_continuous(expression(m[t]), breaks = brk, expand = c(0, 0)) +
  scale_y_continuous("CpG counts", breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.05, 1.12), ylim = c(0, ymax_c), clip = "off") + theme_c + theme(aspect.ratio = 0.95, plot.margin = margin(22, 8, 4, 4)) +
  labs(title = paste0(opt$crukid, "-", s1)) + theme(plot.title = element_text(colour = "grey35", size = 10, hjust = 0.12, margin = margin(b = 14))) +
  ## axis break marker on the truncated first bar
  annotate("segment", x = c(-0.07, -0.07), xend = c(-0.03, -0.03), y = c(ymax_c - 3, ymax_c + 2),
           yend = c(ymax_c + 3, ymax_c + 8), linewidth = 0.5) +
  ann(arr(0.13, 88, 0.095, 64), arr(0.46, 78, 0.44, 58), arr(0.70, 88, 0.84, 44),
      txt(0.46, 86, "?", size = 5))

## ---- 4D ---------------------------------------------------------------------------------
p4D <- ggplot(cd, aes(m_t_1, m_t_2)) + pt() +
  scale_x_continuous(bquote(m[t*","*1]*.(s1)), breaks = brk, expand = c(0, 0)) +
  scale_y_continuous(bquote(m[t*","*2]~.(s1)), breaks = brk, expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.05, 1.1), ylim = c(-0.05, 1.1), clip = "off") + theme_c + theme(aspect.ratio = 1, plot.margin = margin(40, 8, 4, 4)) +
  ann(ell(0.04, 0.77, 0.075, 0.31, cols["mag"]), ell(0.87, 0.85, 0.2, 0.19, cols["grn"]),
      ell(0.69, 0.03, 0.38, 0.09, cols["pur"]),
      txt(0.06, 1.26, "clonal\nmono-allelic", cols["mag"], 3.6, fontface = "bold"),
      txt(0.88, 1.26, "clonal\nbi-allelic", cols["grn"], 3.6, fontface = "bold"),
      txt(0.80, 0.27, "clonal\nmono-allelic", cols["pur"], 3.6, fontface = "bold"))

## ---- 4E ---------------------------------------------------------------------------------
lab_mt <- function(s) bquote(m[t]~.(s))
lim_e <- c(-0.02, 1.08)
e_panel <- function(sx, sy) ggplot(ef, aes(.data[[paste0("m_t_", sx)]], .data[[paste0("m_t_", sy)]])) + pt() +
  scale_x_continuous(lab_mt(sx), breaks = brk, expand = c(0, 0)) +
  scale_y_continuous(lab_mt(sy), breaks = brk, expand = c(0, 0)) +
  coord_cartesian(xlim = lim_e, ylim = lim_e, clip = "off") + theme_g
sB <- samples[2]; sC <- samples[3]                                      # R2, R3
p4E_left <- e_panel(s1, sC) + ggtitle(opt$crukid) +
  theme(plot.title = element_text(colour = "grey35", size = 10, hjust = 0.05)) +
  ann(ell(0.87, 0.86, 0.13, 0.12, cols["grn"]),
      txt(0.2, 0.93, "clonal\nbi-allelic\nDMP", size = 2.9),
      arr(0.36, 0.935, 0.72, 0.89),
      geom_path(data = data.frame(x = c(0.12, 0.40, 0.72, 0.12), y = c(0.06, 0.06, 0.58, 0.06)),
                aes(x, y), colour = "black", linewidth = 0.4, inherit.aes = FALSE),
      txt(0.62, 0.22, "shared clonal\nAS-DMPs and\nsubclonal AS- and\nbi-DMPs", size = 2.9))
p4E_right <- e_panel(s1, sB) +
  ann(ell(0.88, 0.86, 0.14, 0.12, cols["grn"]),
      txt(0.85, 0.24, "distinct\nsubclonal\nAS- and\nbi-DMPs", size = 2.9),
      arr(0.66, 0.20, 0.06, 0.37), arr(0.66, 0.15, 0.40, 0.12))

## ---- 4F ---------------------------------------------------------------------------------
lab_f <- function(s) bquote(m[t*","*1]-m[t*","*2]~.(s))
lim_f <- c(-1.08, 1.08); brk_f <- seq(-1, 1, 0.2)
f_panel <- function(sx, sy) ggplot(ef, aes(.data[[paste0("f_", sx)]], .data[[paste0("f_", sy)]])) + pt() +
  scale_x_continuous(lab_f(sx), breaks = brk_f, expand = c(0, 0)) +
  scale_y_continuous(lab_f(sy), breaks = brk_f, expand = c(0, 0)) +
  coord_cartesian(xlim = lim_f, ylim = lim_f, clip = "off") + theme_g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p4F_left <- f_panel(s1, sC) + ggtitle(opt$crukid) +
  theme(plot.title = element_text(colour = "grey35", size = 10, hjust = 0.05)) +
  ann(ell(0.80, 0.82, 0.17, col = cols["pur"]), ell(0.46, 0.50, 0.17, col = cols["lav"]),
      ell(-0.50, -0.46, 0.17, col = cols["pink"]), ell(-0.80, -0.83, 0.17, col = cols["mag"]),
      txt(-0.55, 0.93, "clonal\nAS-DMP", size = 2.9),
      arr(-0.30, 0.88, 0.62, 0.85), arr(-0.78, 0.75, -0.84, -0.64),
      txt(0.62, -0.52, "shared\nsubclonal\nAS-DMP", size = 2.9),
      arr(0.50, -0.30, 0.46, 0.33), arr(0.30, -0.60, -0.34, -0.50))
p4F_right <- f_panel(s1, sB) +
  ann(ell(0.78, 0.80, 0.17, col = cols["pur"]), ell(0.0, 0.60, 0.17, col = cols["lav"]),
      ell(0.76, 0.0, 0.17, col = cols["lav"]), ell(-0.78, 0.0, 0.17, col = cols["pink"]),
      ell(0.0, -0.62, 0.17, col = cols["pink"]), ell(-0.80, -0.82, 0.17, col = cols["mag"]),
      txt(0.62, -0.70, "distinct\nsubclonal\nAS-DMP", size = 2.9), arr(0.40, -0.80, 0.18, -0.70))

## ---- write ----------------------------------------------------------------------------------
save_both <- function(p, name, w, h) {
  f <- file.path(out_dir, paste0(name, "_", opt$crukid, ".pdf"))
  ggsave(f, p, width = w, height = h)
  ggsave(sub("\\.pdf$", ".png", f), p, width = w, height = h, dpi = 200)
  message("Wrote ", f)
}
p4E <- plot_grid(p4E_left, p4E_right, nrow = 1, align = "h")
p4F <- plot_grid(p4F_left, p4F_right, nrow = 1, align = "h")
save_both(p4C, "Fig4C_mt_histogram", 3.4, 3.3)
save_both(p4D, "Fig4D_mt1_vs_mt2", 3.6, 3.7)
save_both(p4E, "Fig4E_mt_between_regions", 7.2, 3.5)
save_both(p4F, "Fig4F_mdiff_between_regions", 7.2, 3.6)

top <- plot_grid(p4C, p4D, nrow = 1, labels = c("C", "D"), label_size = 18, rel_widths = c(1, 1.05))
fig <- plot_grid(top, p4E, p4F, ncol = 1, labels = c("", "E", "F"), label_size = 18, rel_heights = c(1.05, 1, 1.05))
save_both(fig, "Fig4CDEF", 7.2, 10.8)
