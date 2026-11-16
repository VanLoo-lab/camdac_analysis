#!/usr/bin/env Rscript
#
# Figure 5A: SNV-phased DMRs classified relative to the SNV (mutant) allele.
#   trans      DMR on the WT allele only          cis       DMR on the mutant allele only
#   loss WT    DMR where the WT allele is lost    bi-allelic  DMR on both alleles
#   none       phased bin without a DMR
# each split into hypo / hyper. Left: legend + cartoons (no code existed for these; drawn here).
# Right: pie of all phased bins (n = 4,406) and a zoom of the bins with a DMR (n = 610).
#
# PLOT ONLY. Re-implements, unchanged, the `flag` classification and the two pie plots of
#   Notes_phase_methylation_SNV_v2.R   (scripts/; L885-1146: SNV_DMRs_v2.RData -> flag -> pies)
# and reads the cached per-bin table that script wrote (read-only):
#   ~/rdeocd/txMeth/SNV_DMRs_v2.RData   (`dt`, output of the SNV-phasing / DMR calling at L830-881)
# Nothing is overwritten. The pie rotation (`start`) and the legend/cartoon layout were set by
# hand in the published figure; they are re-created here to match it.
#
#   Rscript plot_5A_SNV_phased_DMR_pies.R [--snv_dmrs_rdata ...] [--out_dir ...]

option_list <- list(
  optparse::make_option("--snv_dmrs_rdata", type = "character", default = "~/rdeocd/txMeth/SNV_DMRs_v2.RData"),
  optparse::make_option("--out_dir", type = "character", default = "~/data/CAMDAC_paper/figure5/outputs"),
  optparse::make_option("--n_AS_DMPs", type = "integer", default = 3L,
                        help = "min. consecutive (AS-)DMPs to call a DMR (Notes: 3)")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(cowplot) })
out_dir <- path.expand(opt$out_dir); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- classification + colours (shared with the 5B-D script; Notes L885-1021, L1126-1131) --------
.script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.script_dir, "fig5_common.R"))
load(path.expand(opt$snv_dmrs_rdata))                                  # dt
dt <- classify_snv_dmrs(dt, opt$n_AS_DMPs)

lv_all <- c("none", "both_hyper", "both_hypo", "LOH_hyper", "LOH_hypo",
            "cis_hyper", "cis_hypo", "trans_hyper", "trans_hypo")      # Notes L1124 (stack order)
tab <- dt[, .N, by = flag][match(lv_all, flag)][, flag := lv_all][is.na(N), N := 0L]
fwrite(tab, file.path(out_dir, "Fig5A_SNV_DMR_phasing_counts.tsv"), sep = "\t")
n_all <- nrow(dt); n_dmr <- dt[flag != "none", .N]
message("phased bins n = ", n_all, "; with DMR n = ", n_dmr)
print(tab)

cl <- snv_dmr_cols
alpha_pie <- 1

pie <- function(d, start, outline = FALSE) {
  d <- d[, .(N = sum(N)), by = flag]
  d[, flag := factor(flag, levels = lv_all[lv_all %in% flag], ordered = TRUE)]
  p <- ggplot(d, aes(x = 1, y = N, fill = flag)) +
    geom_col(width = 2, colour = NA) +
    scale_fill_manual(values = cl, guide = "none") +
    scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_polar("y", start = start) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  if (outline)   # circle outline of the whole population (the "none" slice is ivory)
    p <- p + geom_col(data = data.frame(N = sum(d$N)), aes(x = 1, y = N), width = 2, fill = NA,
                      colour = "black", linewidth = 0.9, inherit.aes = FALSE)
  p
}
# rotations as in the published panel: the DMR wedge is centred on 6 o'clock in the full pie,
# the zoom starts at 3 o'clock (trans-hypo first, bi-allelic hyper last, clockwise).
start_all <- pi - (n_dmr / n_all) * pi                                  # wedge centred at 6 o'clock
p_all <- pie(tab, start_all, outline = TRUE)
p_zoom <- pie(tab[flag != "none"], pi / 2)

## ---- legend + cartoons, drawn in "pixel" coordinates of the published panel (800 x 900) ----------
W <- 800; H <- 828
Y <- function(y) H - y
rrect <- function(x0, x1, y0, y1, r = 9, n = 10) {                       # rounded rectangle polygon
  a <- seq(0, pi / 2, length.out = n)
  cx <- c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a))
  cy <- c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a))
  data.frame(x = cx, y = cy)
}
col_wt <- "#8A2BE2"; fill_wt <- "#C9A9F0"                               # WT chromosome (purple)
col_mut <- "#C8C83C"; fill_mut <- "#F2F2A8"                             # mutant chromosome (yellow)
star_col <- "firebrick"
poly_id <- 0
polys <- list(); pts <- list()
add_poly <- function(d, fill, col) {
  poly_id <<- poly_id + 1; d$id <- poly_id; d$fill <- fill; d$col <- col; polys[[poly_id]] <<- d
}
# one chromosome pair (two segments) at row centre yc (px, y down), left edge x0
chrom <- function(x0, yc, kind) {
  fill <- if (kind == "wt") fill_wt else fill_mut; col <- if (kind == "wt") col_wt else col_mut
  h <- 11
  add_poly(rrect(x0, x0 + 46, Y(yc) - h, Y(yc) + h), fill, col)
  add_poly(rrect(x0 + 50, x0 + 98, Y(yc) - h, Y(yc) + h), fill, col)
}
mark <- function(x, yc, shape) pts[[length(pts) + 1]] <<- data.frame(x = x, y = Y(yc), shape = shape)
cartoon <- function(yc, wt = TRUE, mut = TRUE, wt_dmr = FALSE, mut_dmr = FALSE, snv = TRUE) {
  x0 <- 285
  ywt <- if (mut) yc - 15 else yc; ymu <- if (wt) yc + 15 else yc
  if (wt)  { chrom(x0, ywt, "wt");  if (wt_dmr) mark(x0 + 82, ywt - 6, "dmr") }
  if (mut) { chrom(x0, ymu, "mut"); if (snv) mark(x0 + 66, ymu, "snv"); if (mut_dmr) mark(x0 + 84, ymu + 2, "dmr") }
}
rows_y <- c(trans = 229, cis = 312, loss = 395, bi = 478, none = 598)
cartoon(rows_y["trans"], wt = TRUE, mut = TRUE, wt_dmr = TRUE,  mut_dmr = FALSE)
cartoon(rows_y["cis"],   wt = TRUE, mut = TRUE, wt_dmr = FALSE, mut_dmr = TRUE)
cartoon(rows_y["loss"],  wt = FALSE, mut = TRUE, mut_dmr = TRUE)
cartoon(rows_y["bi"],    wt = TRUE, mut = TRUE, wt_dmr = TRUE,  mut_dmr = TRUE)
cartoon(rows_y["none"],  wt = TRUE, mut = TRUE)
polys_df <- rbindlist(polys); pts_df <- rbindlist(pts)

sq <- function(x, y, s = 14) data.frame(x0 = x - s, x1 = x + s, y0 = Y(y) - s, y1 = Y(y) + s)
sq_df <- rbindlist(list(
  cbind(sq(190, rows_y["trans"]), fill = cl["trans_hypo"]), cbind(sq(240, rows_y["trans"]), fill = cl["trans_hyper"]),
  cbind(sq(190, rows_y["cis"]),   fill = cl["cis_hypo"]),   cbind(sq(240, rows_y["cis"]),   fill = cl["cis_hyper"]),
  cbind(sq(190, rows_y["loss"]),  fill = cl["LOH_hypo"]),   cbind(sq(240, rows_y["loss"]),  fill = cl["LOH_hyper"]),
  cbind(sq(190, rows_y["bi"]),    fill = cl["both_hypo"]),  cbind(sq(240, rows_y["bi"]),    fill = cl["both_hyper"])))
sq_df[, alpha := alpha_pie]

ts <- function(x, y, label, size = 7.6, hjust = 0.5, face = "plain", parse = FALSE, col = "black", angle = 0)
  annotate("text", x = x, y = Y(y), label = label, size = size, hjust = hjust, fontface = face,
           parse = parse, colour = col, angle = angle)

# pie geometry (px) and connectors from the DMR wedge to the zoom pie
c_all <- c(595, 158); r_all <- 132; c_zoom <- c(600, 455); r_zoom <- 125
ang <- function(frac) start_all + 2 * pi * frac                          # clockwise from 12 o'clock
pt_on <- function(cx, cy, r, th) c(cx + r * sin(th), Y(cy) + r * cos(th))
w0 <- pt_on(c_all[1], c_all[2], r_all, ang(0)); w1 <- pt_on(c_all[1], c_all[2], r_all, ang(n_dmr / n_all))
z0 <- pt_on(c_zoom[1], c_zoom[2], r_zoom, -50 * pi / 180)   # connector end points on the zoom pie
z1 <- pt_on(c_zoom[1], c_zoom[2], r_zoom, 52 * pi / 180)

base <- ggplot() +
  geom_polygon(data = polys_df, aes(x, y, group = id, fill = fill, colour = col), linewidth = 0.9) +
  scale_fill_identity() + scale_colour_identity() +
  geom_rect(data = sq_df, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = fill, alpha = alpha),
            colour = NA) + scale_alpha_identity() +
  # hollow square for "none"
  annotate("rect", xmin = 201, xmax = 229, ymin = Y(rows_y["none"]) - 14, ymax = Y(rows_y["none"]) + 14,
           fill = "ivory", colour = "black", linewidth = 1.1) +
  geom_point(data = pts_df[shape == "snv"], aes(x, y), shape = 8, colour = star_col, size = 4.2, stroke = 1.6) +
  geom_point(data = pts_df[shape == "dmr"], aes(x, y), shape = 17, colour = "black", size = 3.6) +
  # header
  ts(160, 118, "hypo", 7.4, hjust = 0, angle = 35) + ts(220, 108, "hyper", 7.4, hjust = 0, angle = 35) +
  annotate("segment", x = 165, xend = 270, y = Y(172), yend = Y(172), linewidth = 0.9) +
  annotate("segment", x = 165, xend = 270, y = Y(545), yend = Y(545), linewidth = 0.9) +
  # row labels
  ts(160, rows_y["trans"], "trans", hjust = 1, face = "italic", size = 8.2) +
  ts(160, rows_y["cis"], "cis", hjust = 1, face = "italic", size = 8.2) +
  ts(160, rows_y["loss"], "loss WT", hjust = 1, size = 8.2) +
  ts(160, rows_y["bi"], "bi-allelic", hjust = 1, size = 8.2) +
  ts(160, rows_y["none"], "none", hjust = 1, size = 8.2) +
  # pie labels + connectors
  annotate("segment", x = w1[1], xend = z0[1], y = w1[2], yend = z0[2], linewidth = 0.7) +
  annotate("segment", x = w0[1], xend = z1[1], y = w0[2], yend = z1[2], linewidth = 0.7) +
  # footer: cartoon key
  ts(48, 688, "cartoon", 6.6, hjust = 0) +
  geom_point(data = data.frame(x = 195, y = Y(688)), aes(x, y), shape = 8, colour = star_col, size = 4.2, stroke = 1.6) +
  ts(212, 688, "SNV", 6.6, hjust = 0) +
  geom_point(data = data.frame(x = 322, y = Y(692)), aes(x, y), shape = 17, size = 3.6) +
  ts(340, 688, "DMR", 6.6, hjust = 0) +
  geom_polygon(data = rbind(cbind(rrect(440, 462, Y(688) - 9, Y(688) + 9, 8), id = 1),
                            cbind(rrect(464, 486, Y(688) - 9, Y(688) + 9, 8), id = 2)),
               aes(x, y, group = id), fill = fill_wt, colour = col_wt, linewidth = 0.9) +
  ts(496, 688, "WT", 6.6, hjust = 0) +
  geom_polygon(data = rbind(cbind(rrect(588, 610, Y(688) - 9, Y(688) + 9, 8), id = 1),
                            cbind(rrect(612, 634, Y(688) - 9, Y(688) + 9, 8), id = 2)),
               aes(x, y, group = id), fill = fill_mut, colour = col_mut, linewidth = 0.9) +
  ts(644, 688, "mutant", 6.6, hjust = 0) +
  # footer: CAMDAC methylation key
  ts(400, 748, "CAMDAC methylation", 6.9) +
  geom_point(data = data.frame(x = c(95, 195, 515, 688), y = Y(803), col = c("black", "rosybrown", "#CDCD42", "#7B2FD6")),
             aes(x, y, colour = col), size = 5.5, shape = 16) +
  geom_point(data = data.frame(x = 300, y = Y(803)), aes(x, y), shape = 15, colour = "grey60", size = 5) +
  ts(112, 803, "m[t]", 7.4, hjust = 0, parse = TRUE) + ts(212, 803, "m[n]", 7.4, hjust = 0, parse = TRUE) +
  ts(315, 803, "99% HDI", 6.6, hjust = 0) + ts(530, 803, "mutant", 6.6, hjust = 0) + ts(703, 803, "WT", 6.6, hjust = 0) +
  ts(18, 42, "A", 11, hjust = 0, face = "bold") +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) + theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# coord_polar draws the circle with diameter 0.8 of the panel, so the half-width of the panel is r / 0.8
place <- function(p, cx, cy, r) { R <- r / 0.8
  draw_plot(p, x = (cx - R) / W, y = (H - (cy + R)) / H, width = 2 * R / W, height = 2 * R / H) }
lab <- function(x, y, txt) draw_label(txt, x = x / W, y = (H - y) / H, size = 24, fontface = "plain")
fig <- ggdraw(base) + place(p_all, c_all[1], c_all[2], r_all) + place(p_zoom, c_zoom[1], c_zoom[2], r_zoom) +
  lab(c_all[1] - 5, 95, paste0("n= ", format(n_all, big.mark = ","))) + lab(c_zoom[1] - 5, 628, paste0("n= ", n_dmr))

ggsave(file.path(out_dir, "Fig5A_SNV_phased_DMR_pies.pdf"), fig, width = 8, height = 8.28)
ggsave(file.path(out_dir, "Fig5A_SNV_phased_DMR_pies.png"), fig, width = 8, height = 8.28, dpi = 150)
message("Wrote ", file.path(out_dir, "Fig5A_SNV_phased_DMR_pies.pdf"))
