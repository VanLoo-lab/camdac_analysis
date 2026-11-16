#!/usr/bin/env Rscript
#
# Figure 5B-D: example SNV-phased DMRs (one bin of a tumour region, CpG by CpG).
#   B  CDKN2A, CRUK0082-R2 (CRUK0082)  cis   hyper (DMR on the mutant allele)
#   C  MDGA1,  CRUK0071-R7 (CRUK0071)  trans hyper (DMR on the WT allele)
#   D  CDH10,  CRUK0072-R4 (CRUK0072)  loss-of-WT hypo (WT allele lost)
#
# Per CpG: black = purified tumour methylation m_t, rosybrown = normal m_n (points and 99% HDI ribbon;
# grey ribbon = m_t HDI), purple / yellow = SNV-phased m_t,1 (WT allele) / m_t,2 (mutant allele),
# red star = SNV. The black "composite" dots at phased CpGs are (m_t,1*CN_1 + m_t,2*CN_2)/(CN_1 + CN_2).
#
# PLOT ONLY. Plot code is the `Methylation_plot` block of Notes_phase_methylation_SNV_v2.R
# (scripts/, L1296-1450; the 2022 variant in R_scripts/Notes_phase_methylation_SNV_1_plus_0.R
# L2031-2146 is the same plot with circles), with the published styling. Reads (read-only):
#   ~/rdeocd/txMeth/SNV_DMRs_v2.RData                                 DMR bins, to class + bin coordinates
#   ~/rdeocd/pipeline_files/hg19_annotations/hg19_all_regions_annotations.RData   cluster_id -> start/end
#   ~/rdeocd/txMeth/<patient>/Methylation/<R>/re_segmented/purified_tumour_DMP_calls.RData
#   ~/rdeocd/txMeth/<patient>/SNVs/<R>/SNV_phased_DMRs_calls_<R>.RData
# The class (square next to the gene name and the allele cartoon) is derived with the same
# classification as Fig 5A (fig5_common.R). Add other loci with --examples "gene:patient:region[:class],...".
# The optional 4th field overrides the derived class: the Aug-2021 cache has m_t_2_DMPs = NA for the
# CDH10 bin (LOH, CN_1 = 0), so the classification returns "none" there although the bin is a
# loss-of-WT hypo DMR (as in the paper), hence `LOH_hypo` is given explicitly for panel D.

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
  optparse::make_option("--snv_dmrs_rdata", type = "character", default = "~/rdeocd/txMeth/SNV_DMRs_v2.RData"),
  optparse::make_option("--annotations", type = "character",
    default = "~/rdeocd/pipeline_files/hg19_annotations/hg19_all_regions_annotations.RData"),
  optparse::make_option("--out_dir", type = "character", default = "~/data/CAMDAC_paper/figure5/outputs"),
  optparse::make_option("--examples", type = "character",
    default = "CDKN2A:CRUK0082:R2,MDGA1:CRUK0071:R7,CDH10:CRUK0072:R4:LOH_hypo")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.script_dir, "fig5_common.R"))
tx_dir <- path.expand(opt$tx_dir)
out_dir <- path.expand(opt$out_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- bins (class + coordinates) -----------------------------------------------------------
load(path.expand(opt$snv_dmrs_rdata))                                          # dt
dt <- classify_snv_dmrs(dt)
e <- new.env(); load(path.expand(opt$annotations), envir = e)
cl_xy <- unique(e$all_regions_anno[, .(cluster_id, chr = as.character(chrom), start, end)])
rm(e)
dt <- merge(dt, cl_xy, by = "cluster_id", all.x = TRUE, suffixes = c("", ".cl"))
clin <- fread(file.path(tx_dir, "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))

## ---- cartoon / gene-label helpers (pixel coordinates) -----------------------------------------
rrect <- function(x0, x1, y0, y1, r = 5, n = 8) {
  a <- seq(0, pi / 2, length.out = n)
  data.frame(x = c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a)),
             y = c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a)))
}
col_wt <- "#8A2BE2"; fill_wt <- "#C9A9F0"; col_mut <- "#C8C83C"; fill_mut <- "#F2F2A8"
# allele cartoon for a class: trans = DMR on WT; cis = DMR on mutant; LOH = only the mutant allele; both = both
cartoon_gg <- function(cls) {
  kind <- sub("_.*", "", cls)
  has_wt <- kind != "LOH"; wt_dmr <- kind %in% c("trans", "both"); mut_dmr <- kind %in% c("cis", "LOH", "both")
  polys <- list(); pts <- list()
  chrom <- function(yc, wt) {
    f <- if (wt) fill_wt else fill_mut; cc <- if (wt) col_wt else col_mut
    for (x0 in c(2, 38)) polys[[length(polys) + 1]] <<- cbind(rrect(x0, x0 + 33, yc - 8, yc + 8), id = length(polys) + 1, f = f, c = cc)
  }
  ywt <- 44; ymu <- if (has_wt) 22 else 33
  if (has_wt) { chrom(ywt, TRUE); if (wt_dmr) pts[[1]] <- data.frame(x = 64, y = ywt + 7, s = 17) }
  chrom(ymu, FALSE)
  pts[[length(pts) + 1]] <- data.frame(x = 50, y = ymu, s = 8)
  if (mut_dmr) pts[[length(pts) + 1]] <- data.frame(x = 62, y = ymu + 1, s = 17)
  P <- rbindlist(polys); S <- rbindlist(pts)
  ggplot() + geom_polygon(data = P, aes(x, y, group = id, fill = f, colour = c), linewidth = 0.7) +
    scale_fill_identity() + scale_colour_identity() +
    geom_point(data = S[s == 8], aes(x, y), shape = 8, colour = "firebrick", size = 2.3, stroke = 0.9) +
    geom_point(data = S[s == 17], aes(x, y), shape = 17, colour = "black", size = 2) +
    coord_fixed(xlim = c(0, 76), ylim = c(8, 58), expand = FALSE) + theme_void()
}
gene_label_gg <- function(gene, cls) {
  w <- 28 + 9 * nchar(gene)
  ggplot() + geom_polygon(data = rrect(2, w + 36, 4, 36, 8), aes(x, y), fill = "white", colour = "black", linewidth = 0.6) +
    annotate("text", x = 10, y = 20, label = gene, hjust = 0, size = 3.7) +
    annotate("rect", xmin = w + 8, xmax = w + 30, ymin = 9, ymax = 31, fill = snv_dmr_cols[[cls]]) +
    coord_fixed(xlim = c(0, w + 38), ylim = c(0, 40), expand = FALSE) + theme_void()
}

## ---- one panel ---------------------------------------------------------------------------------
make_panel <- function(gene, patient, region, letter, cls_override = NA) {
  # `patient` may be CRUKID; resolve to internal Patient for paths / LTX_ID
  cruk_in <- patient
  patient <- resolve_tx_patient_id(patient, .tx_map)
  crukid <- clin[Patient == patient, CRUKID][1]
  if (is.na(crukid) || !nzchar(crukid)) crukid <- cruk_in
  lab <- paste0(sub("^CRUK0*", "", crukid), "-", region)
  pat_re <- paste0("^", gene, "|;", gene)
  bin <- dt[LTX_ID == paste(patient, region, sep = "_") & (grepl(pat_re, gene_ids) | grepl(pat_re, gene_prom_ids))][1]
  stopifnot(nrow(bin) == 1)
  cls <- bin$flag
  if (!is.na(cls_override)) {
    if (cls != cls_override) message("   derived class '", cls, "' replaced by '", cls_override, "'")
    cls <- cls_override
  }
  message(letter, ": ", gene, " ", lab, " bin chr", bin$chr, ":", bin$start, "-", bin$end, "  class ", cls)

  em <- new.env()
  load(file.path(tx_dir, patient, "Methylation", region, "re_segmented", "purified_tumour_DMP_calls.RData"), envir = em)
  df_anno <- as.data.table(em$df_mt)[as.character(chromosome) == bin$chr & start >= bin$start & end <= bin$end]
  ep <- new.env()
  load(file.path(tx_dir, patient, "SNVs", region, paste0("SNV_phased_DMRs_calls_", region, ".RData")), envir = ep)
  df_ph <- as.data.table(ep$df_phased_methylation)[as.character(chrom) == bin$chr & start >= bin$start & end <= bin$end]
  snv_pos <- unique(as.numeric(df_ph$POS_genotype))
  cap <- function(x) ifelse(is.na(x), NA_real_, pmin(pmax(x, 0), 1))
  message("   ", nrow(df_anno), " CpGs, ", nrow(df_ph), " phased, SNV at ", paste(snv_pos, collapse = ","))

  kb <- (bin$end - bin$start) > 300                       # kb axis for the longer bins, bp for the short one (as published)
  xl <- if (kb) "Chromosome position (in kb)" else "Chromosome position (in bp)"
  xs <- if (kb) scale_x_continuous(breaks = scales::breaks_width(200), labels = function(x) format(x / 1000, nsmall = 1),
                                   expand = expansion(mult = 0.03))
        else scale_x_continuous(breaks = scales::breaks_width(25), labels = function(x) format(x, scientific = FALSE),
                                expand = expansion(mult = 0.05))
  p <- ggplot(df_anno) + theme_minimal(base_size = 11) +
    geom_ribbon(aes(x = start, ymin = cap(m_t_low), ymax = cap(m_t_high)), alpha = 0.5, fill = "lightgrey") +
    geom_ribbon(aes(x = start, ymin = m_n_low, ymax = m_n_high), alpha = 0.5, fill = "rosybrown") +
    geom_point(aes(x = start, y = m_n), size = 2.2, colour = "rosybrown", shape = 16) +
    geom_point(aes(x = start, y = m_t_corr), size = 2.2, colour = "black", shape = 16)
  if (nrow(df_ph) > 0)
    p <- p +
      geom_point(data = df_ph, aes(x = start, y = cap(m_t_1), colour = "WT"), size = 2.2, shape = 16, na.rm = TRUE) +
      geom_point(data = df_ph, aes(x = start, y = cap(m_t_2), colour = "mutant"), size = 2.2, shape = 16) +
      geom_point(data = df_ph, aes(x = start, y = (m_t_1 * CN_1 + m_t_2 * CN_2) / (CN_1 + CN_2)),
                 size = 2.2, colour = "black", shape = 16, na.rm = TRUE)
  p <- p +
    scale_colour_manual(values = c(mutant = "#CDCD00", WT = "#7B2FD6"), guide = "none") +
    annotate("point", x = snv_pos, y = -0.1, shape = 8, colour = "firebrick", size = 4, stroke = 1.3) +
    xs + scale_y_continuous("Methylation rate", breaks = seq(0, 1, 0.2), labels = function(x) format(x, nsmall = 1)) +
    coord_cartesian(ylim = c(-0.13, 1.06), clip = "off") + xlab(xl) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.ticks = element_line(colour = "black"),
          panel.grid.major = element_line(colour = "grey88", linewidth = 0.6), panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black", size = 9), axis.title = element_text(size = 11),
          plot.margin = margin(2, 16, 2, 4))
  # header: sample label above, gene (+class square) and allele cartoon
  header <- ggdraw() + draw_label(lab, x = 0.5, y = 0.84, size = 11) +
    draw_plot(gene_label_gg(gene, cls), x = 0.20, y = 0.0, width = 0.45, height = 0.7, hjust = 0, vjust = 0) +
    draw_plot(cartoon_gg(cls), x = 0.70, y = 0.0, width = 0.28, height = 0.7)
  list(fig = plot_grid(header, p, ncol = 1, rel_heights = c(0.2, 1)), cls = cls, bin = bin)
}

ex <- strsplit(strsplit(opt$examples, ",")[[1]], ":")
letters_ <- c("B", "C", "D", "E", "F", "G")[seq_along(ex)]
panels <- Map(function(x, l) make_panel(x[1], x[2], x[3], l, if (length(x) > 3) x[4] else NA), ex, letters_)
hts <- c(1, 1.08, 1)[seq_along(ex)]
for (i in seq_along(ex)) {
  f <- file.path(out_dir, sprintf("Fig5%s_SNV_phased_DMR_%s_%s-%s", letters_[i], ex[[i]][1], ex[[i]][2], ex[[i]][3]))
  g <- plot_grid(panels[[i]]$fig, labels = letters_[i], label_size = 18)
  ggsave(paste0(f, ".pdf"), g, width = 3.6, height = 3.3); ggsave(paste0(f, ".png"), g, width = 3.6, height = 3.3, dpi = 200)
}
all_fig <- plot_grid(plotlist = lapply(panels, `[[`, "fig"), ncol = 1, labels = letters_, label_size = 18,
                     rel_heights = hts)
ggsave(file.path(out_dir, "Fig5BCD_SNV_phased_DMR_examples.pdf"), all_fig, width = 3.6, height = 3.3 * length(ex))
ggsave(file.path(out_dir, "Fig5BCD_SNV_phased_DMR_examples.png"), all_fig, width = 3.6, height = 3.3 * length(ex), dpi = 200)
message("Wrote ", file.path(out_dir, "Fig5BCD_SNV_phased_DMR_examples.pdf"))
