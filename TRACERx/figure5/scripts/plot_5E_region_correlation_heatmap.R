#!/usr/bin/env Rscript
#
# Figure 5E: pairwise Pearson correlations of methylation between the tumour regions of one patient
# (CRUK0062, regions R1-R7) and the matched normal, at the patient's DMP loci.
#   upper-left triangle  m_b  bulk methylation (loci that are bulk DMPs in >= 1 region)
#   lower-right triangle m_t  CAMDAC-purified methylation (loci that are purified DMPs in >= 1 region)
#   top row / right column: correlation of each region with the normal (m_b / m_t, normal: m_b = m_t = m_n)
# Regions are ordered by tumour purity (rho_reseg), highest purity at the bottom / left.
#
# Method = Notes_clustering_DMPs_MVPs_pilot.R L19-165 (the `corr_mb` / `corr_mt` loop and the two
# scale_fill_gradient2() heatmaps; the 2019 per-patient output was Clustering/<Patient>_m{b,t}_corr_all_DMPs.pdf):
#   m_diff_t = m_t_corr - m_n ; DMP_t = prob_hyper_corr < 0.01 & m_diff_t >= 0.2 (hyper) | prob_hypo_corr < 0.01 &
#   m_diff_t <= -0.2 (hypo);  DMP_b the same with m_b ; pairwise.complete Pearson correlation.
# The split-cell layout, colour bars and axes are re-drawn here to match the published panel.
#
# Input (read-only): by default the per-region tables saved by the published analysis
# (~/rdeocd/txMeth/<patient>/Clustering/clustering_workspace_<patient>.RData, `list_patient`), which reproduces
# its df_results exactly. The current Methylation/<R>/purified_tumour_DMP_calls.RData files have been re-run since
# (different m_t_corr / DMP sets, m_t-vs-normal correlations ~0.3 higher), use --clustering_workspace none for those.
# Writes
#   outputs/Fig5E_region_correlations.tsv, Fig5E_region_correlation_heatmap.{pdf,png}
#
#   Rscript plot_5E_region_correlation_heatmap.R [--patient CRUK0062] [--qc_workspace <clustering_workspace.RData>]

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
  optparse::make_option("--patient", type = "character", default = "CRUK0062"),
  optparse::make_option("--out_dir", type = "character", default = "~/data/CAMDAC_paper/figure5/outputs"),
  optparse::make_option("--recompute", type = "logical", default = FALSE,
                        help = "ignore the cached outputs/Fig5E_region_correlations.tsv"),
  optparse::make_option("--clustering_workspace", type = "character", default = "auto",
    help = paste("RData with the saved per-region tables (list_patient) of the published analysis; 'auto' =",
                 "<tx_dir>/<patient>/Clustering/clustering_workspace_<patient>.RData; 'none' = recompute DMP calls from",
                 "the current purified_tumour_DMP_calls.RData files (their m_t_corr values have changed since 2019)")),
  optparse::make_option("--qc_workspace", type = "character", default = "",
    help = "optional: 2019 clustering_workspace_<patient>.RData; prints the max |difference| to its df_results")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (!is.null(opt$patient) && nzchar(opt$patient))
  opt$patient <- resolve_tx_patient_id(opt$patient, .tx_map)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
tx_dir <- path.expand(opt$tx_dir); out_dir <- path.expand(opt$out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tsv <- file.path(out_dir, "Fig5E_region_correlations.tsv")

clin <- fread(file.path(tx_dir, "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))
cl <- clin[Patient == opt$patient & !grepl("LN", Sample)]
cl <- cl[order(-rho_reseg)]                                                  # highest purity first
regions <- cl$Sample
message(opt$patient, " (", cl$CRUKID[1], "): regions by purity: ", paste0(regions, " (", cl$rho_reseg, ")", collapse = ", "))

## ---- compute (pilot script L19-135) ------------------------------------------------------------
if (!file.exists(tsv) || opt$recompute) {
  dat <- list(); dmp_b <- list(); dmp_t <- list()
  ws <- if (opt$clustering_workspace == "auto") file.path(tx_dir, opt$patient, "Clustering", paste0("clustering_workspace_", opt$patient, ".RData")) else opt$clustering_workspace
  use_ws <- !(ws %in% c("none", "")) && file.exists(path.expand(ws))
  if (use_ws) {                                  # per-region tables saved by the analysis behind the published figure
    message("Using the saved per-region tables (list_patient) of ", ws)
    ew <- new.env(); load(path.expand(ws), envir = ew)
    lp <- setNames(ew$list_patient[seq_along(regions)], NULL)
    names(lp) <- vapply(lp, function(x) sub("^CRUK[0-9]+_", "", as.character(x$label[1])), "")
    rm(ew)
  } else message("Recomputing the DMP calls from the current purified_tumour_DMP_calls.RData files")
  for (r in regions) {
    if (use_ws) {
      d <- as.data.table(lp[[r]])[, .(chromosome = as.character(chromosome), start, end, m_b, m_n, m_t, DMP_t, DMP_b)]
    } else {
      e <- new.env()
      load(file.path(tx_dir, opt$patient, "Methylation", r, "purified_tumour_DMP_calls.RData"), envir = e)
      d <- as.data.table(e$df_mt)[, .(chromosome = as.character(chromosome), start, end, m_b, m_n, m_t = m_t_corr,
                                      prob_hyper_corr, prob_hypo_corr)]
      rm(e)
      d[, `:=`(m_diff_t = m_t - m_n, m_diff_b = m_b - m_n)]
      d[, DMP_t := ifelse(prob_hyper_corr < 0.01 & m_diff_t >= 0.2, "hyper",
                   ifelse(prob_hypo_corr < 0.01 & m_diff_t <= -0.2, "hypo", NA_character_))]
      d[, DMP_b := ifelse(prob_hyper_corr < 0.01 & m_diff_b >= 0.2, "hyper",
                   ifelse(prob_hypo_corr < 0.01 & m_diff_b <= -0.2, "hypo", NA_character_))]
    }
    d <- d[!duplicated(d[, .(chromosome, start, end)])]
    dmp_b[[r]] <- d[!is.na(DMP_b), .(chromosome, start, end)]
    dmp_t[[r]] <- d[!is.na(DMP_t), .(chromosome, start, end)]
    dat[[r]] <- d[, .(chromosome, start, end, m_b, m_t, m_n)]
    message(r, ": ", nrow(d), " CpGs, ", nrow(dmp_b[[r]]), " bulk DMPs, ", nrow(dmp_t[[r]]), " purified DMPs")
  }
  # normal profile: m_n at all loci seen in any region (m_b = m_t = m_n)
  nrm <- unique(rbindlist(lapply(dat, function(d) d[, .(chromosome, start, end, m_n)])))
  nrm <- nrm[!duplicated(nrm[, .(chromosome, start, end)])]
  loci_b <- unique(rbindlist(dmp_b)); loci_t <- unique(rbindlist(dmp_t))
  setkey(loci_b, chromosome, start, end); setkey(loci_t, chromosome, start, end)
  getmat <- function(loci, col) {
    M <- loci[, .(chromosome, start, end)]
    for (r in regions) M[[r]] <- dat[[r]][loci, on = .(chromosome, start, end)][[col]]
    M[["Normal"]] <- nrm[loci, on = .(chromosome, start, end)][["m_n"]]
    as.matrix(M[, c(regions, "Normal"), with = FALSE])
  }
  Cb <- cor(getmat(loci_b, "m_b"), use = "pairwise.complete.obs", method = "pearson")
  Ct <- cor(getmat(loci_t, "m_t"), use = "pairwise.complete.obs", method = "pearson")
  nm <- c(regions, "Normal")
  res <- CJ(sample1 = nm, sample2 = nm, sorted = FALSE)
  res[, `:=`(corr_mb = Cb[cbind(sample1, sample2)], corr_mt = Ct[cbind(sample1, sample2)])]
  res[, `:=`(n_loci_mb = nrow(loci_b), n_loci_mt = nrow(loci_t))]
  fwrite(res, tsv, sep = "\t")
}
res <- fread(tsv)

if (nzchar(opt$qc_workspace)) {                                              # QC vs the 2019 workspace
  e <- new.env(); load(path.expand(opt$qc_workspace), envir = e)
  old <- as.data.table(e$df_results); old[, c("sample1", "sample2") := lapply(.SD, function(x) sub("^CRUK[0-9]+_", "", as.character(x))), .SDcols = c("sample1", "sample2")]
  m <- merge(res, old, by = c("sample1", "sample2"), suffixes = c("", "_2019"))
  message("QC vs 2019 workspace: ", nrow(m), " pairs; max |d corr_mb| = ", signif(max(abs(m$corr_mb - m$corr_mb_2019)), 3),
          ", max |d corr_mt| = ", signif(max(abs(m$corr_mt - m$corr_mt_2019)), 3))
}

## ---- drawing (all in "cell" units) -------------------------------------------------------------------
lim_min <- min(c(res$corr_mt, res$corr_mb), na.rm = TRUE)                      # pilot L143-146
lim_min <- seq(-1, 0, 0.05)[which.min(abs(seq(-1, 0, 0.05) - lim_min))]
if (lim_min > min(c(res$corr_mt, res$corr_mb), na.rm = TRUE)) lim_min <- lim_min - 0.05
pal_b <- scales::div_gradient_pal("chartreuse3", "white", "mediumpurple3", "Lab")      # scale_fill_gradient2 (pilot L151)
pal_t <- scales::div_gradient_pal("red", "white", "steelblue", "Lab")                  # (pilot L157)
fill_of <- function(x, pal) pal(scales::rescale_mid(x, to = c(0, 1), from = c(lim_min, 1), mid = 0))

nR <- length(regions)
xs <- setNames(seq_len(nR), regions)                       # columns left -> right: highest purity first
ys <- setNames(seq_len(nR), regions)                       # row y: highest purity at the bottom (y = 1), R6 on top
xN <- nR + 1.25; yN <- nR + 1.25                           # Normal column / row centres
val <- function(a, b, what) res[sample1 == a & sample2 == b, get(what)]
fmt <- function(x) sub("^-0\\.00$", "0.00", sprintf("%.2f", x))

tiles <- list(); tri <- list(); labs <- list()
for (a in regions) for (b in regions) {                    # a = column region (x), b = row region (y)
  x <- xs[[a]]; y <- ys[[b]]
  if (a == b) {                                            # identity cell: split along the anti-diagonal
    tri[[length(tri) + 1]] <- data.frame(x = x + c(-.5, .5, -.5), y = y + c(.5, .5, -.5), fill = fill_of(1, pal_b), id = length(tri) + 1)
    tri[[length(tri) + 1]] <- data.frame(x = x + c(.5, .5, -.5), y = y + c(.5, -.5, -.5), fill = fill_of(1, pal_t), id = length(tri) + 1)
  } else {
    # upper-left half of the matrix (column more pure than row ... i.e. y > x) -> m_b ; lower-right -> m_t
    wh <- if (y > x) "corr_mb" else "corr_mt"
    v <- val(a, b, wh)
    tiles[[length(tiles) + 1]] <- data.frame(x = x, y = y, w = 1, h = 1, fill = fill_of(v, if (wh == "corr_mb") pal_b else pal_t))
    labs[[length(labs) + 1]] <- data.frame(x = x, y = y, lab = fmt(v))
  }
}
# top row: Normal vs region (m_b); right column: region vs Normal (m_t)
for (a in regions) {
  v <- val(a, "Normal", "corr_mb")
  tiles[[length(tiles) + 1]] <- data.frame(x = xs[[a]], y = yN, w = 1, h = 0.8, fill = fill_of(v, pal_b))
  labs[[length(labs) + 1]] <- data.frame(x = xs[[a]], y = yN, lab = fmt(v))
}
for (b in regions) {
  v <- val(b, "Normal", "corr_mt")
  tiles[[length(tiles) + 1]] <- data.frame(x = xN, y = ys[[b]], w = 0.8, h = 1, fill = fill_of(v, pal_t))
  labs[[length(labs) + 1]] <- data.frame(x = xN, y = ys[[b]], lab = fmt(v))
}
tri[[length(tri) + 1]] <- data.frame(x = xN + .4 * c(-1, 1, -1), y = yN + .4 * c(1, 1, -1), fill = fill_of(1, pal_b), id = length(tri) + 1)
tri[[length(tri) + 1]] <- data.frame(x = xN + .4 * c(1, 1, -1), y = yN + .4 * c(1, -1, -1), fill = fill_of(1, pal_t), id = length(tri) + 1)
T_ <- rbindlist(tiles); TR <- rbindlist(tri); L_ <- rbindlist(labs)
# the Normal-vs-lowest-corr cell is white in the paper (0.00): fill_of handles it (colour ~ white)

lw <- 1.6
ax <- function(...) annotate("segment", ..., linewidth = lw, lineend = "square")
tick <- 0.12
p <- ggplot() +
  geom_tile(data = T_, aes(x, y, width = w, height = h, fill = fill)) +
  geom_polygon(data = TR, aes(x, y, group = id, fill = fill)) +
  scale_fill_identity() +
  geom_text(data = L_, aes(x, y, label = lab), size = 4.4) +
  # axes (thick black lines + ticks + labels), as in the figure
  ax(x = 0.4, xend = 0.4, y = 0.5, yend = nR + 0.5) + ax(x = 0.4, xend = 0.4, y = yN - 0.4, yend = yN + 0.4) +
  ax(x = 0.5, xend = nR + 0.5, y = 0.35, yend = 0.35) + ax(x = xN - 0.4, xend = xN + 0.4, y = 0.35, yend = 0.35) +
  annotate("segment", x = 0.4 - tick, xend = 0.4, y = ys, yend = ys, linewidth = lw / 1.5) +
  annotate("segment", x = 0.4 - tick, xend = 0.4, y = yN, yend = yN, linewidth = lw / 1.5) +
  annotate("segment", x = xs, xend = xs, y = 0.35 - tick, yend = 0.35, linewidth = lw / 1.5) +
  annotate("segment", x = xN, xend = xN, y = 0.35 - tick, yend = 0.35, linewidth = lw / 1.5) +
  annotate("text", x = 0.4 - tick - 0.08, y = ys, label = names(ys), hjust = 1, size = 7.2, colour = "grey25") +
  annotate("text", x = 0.4 - tick - 0.08, y = yN, label = "Normal", hjust = 1, size = 7.2, colour = "grey25") +
  annotate("text", x = xs, y = 0.35 - tick - 0.12, label = names(xs), vjust = 1, size = 7.2, colour = "grey25") +
  annotate("text", x = xN, y = 0.35 - tick - 0.12, label = "Normal", vjust = 1, size = 7.2, colour = "grey25") +
  # purity wedge + label
  annotate("polygon", x = c(-1.65, -1.95, -1.35), y = c(nR + 0.4, 0.62, 0.62), fill = "grey45") +
  annotate("text", x = -2.2, y = (nR + 1) / 2, label = "Fraction tumor DNA", angle = 90, size = 7.6) +
  coord_fixed(xlim = c(-2.6, xN + 1.1), ylim = c(-3.1, yN + 0.7), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

## ---- colour bars (bottom) ----------------------------------------------------------------------------
bar <- function(x0, y0, pal, title, tcol = "black", w = 2.7, h = 0.36) {
  v <- seq(1, lim_min, length.out = 120)                           # 1.0 on the left -> lim_min on the right
  d <- data.frame(x = x0 + (seq_along(v) - 0.5) / length(v) * w, y = y0, fill = fill_of(v, pal))
  xpos <- function(val) x0 + (1 - val) / (1 - lim_min) * w
  up <- c(1.0, 0.6, 0.2, -0.2); dn <- c(0.8, 0.4, 0.0)
  list(geom_tile(data = d, aes(x, y, fill = fill), width = w / length(v) + 0.002, height = h),
       annotate("segment", x = xpos(c(up, dn)), xend = xpos(c(up, dn)), y = y0 - h / 2, yend = y0 - h / 2 + 0.08, colour = "white", linewidth = 0.5),
       annotate("segment", x = xpos(c(up, dn)), xend = xpos(c(up, dn)), y = y0 + h / 2, yend = y0 + h / 2 - 0.08, colour = "white", linewidth = 0.5),
       annotate("text", x = xpos(up), y = y0 + h / 2 + 0.2, label = ifelse(up < 0, paste0("\u2212", abs(up)), format(up, nsmall = 1)), size = 5.2),
       annotate("text", x = xpos(dn), y = y0 - h / 2 - 0.2, label = format(dn, nsmall = 1), size = 5.2),
       annotate("text", x = x0 - 0.15, y = y0, label = title, hjust = 1, size = 7.4, fontface = "bold", parse = TRUE))
}
yb <- -2.1
p <- p + bar(1.4, yb, pal_t, "m[t]") + bar(5.7, yb, pal_b, "m[b]") +
  annotate("text", x = -2.55, y = yb + 0.05, label = "Pearson\ncorrelations", hjust = 0, size = 5.6, lineheight = 0.95)

ggsave(file.path(out_dir, "Fig5E_region_correlation_heatmap.pdf"), p, width = 7.2, height = 7.9, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig5E_region_correlation_heatmap.png"), p, width = 7.2, height = 7.9, dpi = 150)
message("Wrote ", file.path(out_dir, "Fig5E_region_correlation_heatmap.pdf"))
print(dcast(res, sample1 ~ sample2, value.var = "corr_mb")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
