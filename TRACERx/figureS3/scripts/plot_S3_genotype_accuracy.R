#' @title Supplementary Figure S3A–D: RRBS genotyping accuracy vs WGS
#' @author Lili Larose Cadieux
#' @description Heterozygous SNP calling from normal RRBS vs WGS ground truth
#' (pilot TRACERx patients). Panels:
#' * A - false-positive rates by SNP type / sequence context + BAF-error hist
#' * B - false-negative rates by SNP type / sequence context + BAF-error hist
#' * C - FNR vs FPR across coverage and heterozygosity boundaries
#' * D - TP/FN enrichment at CCGG-near / short-read / low-mq loci (CRUK0062)
#'
#' Extracted from \code{Notes_RRBS_WGS_Genotype_LogR_BAF_comparison.R}
#' (block marked "paper figures 3 and S3"). Plots from the cached pilot object
#' \code{all_pilot_patients_genotype_TNR_FNR.RData} (no WGS re-comparison).
#'
#' S3E–F (WGBS / CRC) live under
#' \code{camdac_analysis/CRC/Figure S3 Genotyping Accuracy/} and are not
#' generated here.

x <- c("data.table", "ggplot2", "optparse")
invisible(lapply(x, function(y) {
  suppressWarnings(
    suppressMessages(
      library(y, character.only = TRUE, quietly = TRUE,
              warn.conflicts = FALSE)
    )
  )
}))
rm(x)

################################################################################
# Helpers
################################################################################

resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) {
    return(path)
  }
  file.path(base_dir, path)
}

WIDTH_LABELS <- c(
  "1" = "Other SNPs",
  "2" = "SNP at CpG",
  "4" = "SNP at CCGG"
)

collapse_snp_type <- function(ref, alt) {
  snp <- paste0(ref, ">", alt)
  ifelse(
    snp %in% c("A>C", "C>A", "G>T", "T>G"), "C>A",
    ifelse(
      snp %in% c("A>G", "G>A", "C>T", "T>C"), "C>T",
      ifelse(
        snp %in% c("A>T", "T>A"), "T>A",
        ifelse(snp %in% c("C>G", "G>C"), "C>G", snp)
      )
    )
  )
}

reclassify_type <- function(dt, thres = 0.1) {
  d <- copy(dt)
  d[, type := ifelse(
    is.finite(BAF_n) & BAF_n >= thres & BAF_n <= (1 - thres),
    "Heterozygous", "Homozygous"
  )]
  d[]
}

summarise_fpr <- function(df2) {
  # WGS-true homozygous: FP = called heterozygous in RRBS
  d <- as.data.table(df2)
  d[, width := as.character(width)]
  d[, new_SNP := collapse_snp_type(ref, alt)]
  out <- d[
    ,
    .(
      Het_SNPs = sum(type == "Heterozygous", na.rm = TRUE),
      Hom_SNPs = sum(type == "Homozygous", na.rm = TRUE)
    ),
    by = .(new_SNP, width)
  ]
  out <- out[new_SNP %in% c("C>A", "C>G", "C>T", "T>A")]
  out[, FPR := 100 * Het_SNPs / (Het_SNPs + Hom_SNPs)]
  out[, width := factor(width, levels = c("1", "2", "4"), labels = WIDTH_LABELS)]
  out[]
}

summarise_fnr <- function(df1) {
  # WGS-true heterozygous: FN = called homozygous in RRBS
  d <- as.data.table(df1)
  d[, width := as.character(width)]
  d[, new_SNP := collapse_snp_type(ref, alt)]
  out <- d[
    ,
    .(
      Het_SNPs = sum(type == "Heterozygous", na.rm = TRUE),
      Hom_SNPs = sum(type == "Homozygous", na.rm = TRUE)
    ),
    by = .(new_SNP, width)
  ]
  out <- out[new_SNP %in% c("C>A", "C>G", "C>T", "T>A")]
  out[, FNR := 100 * Hom_SNPs / (Het_SNPs + Hom_SNPs)]
  out[, width := factor(width, levels = c("1", "2", "4"), labels = WIDTH_LABELS)]
  out[]
}

plot_fpr_bubbles <- function(df2, outfile) {
  df_summary <- summarise_fpr(df2)
  p <- ggplot(df_summary, aes(x = new_SNP, y = 1)) +
    geom_point(aes(size = Het_SNPs, colour = FPR)) +
    geom_text(aes(label = round(FPR, 1)), size = 3.5, colour = "black") +
    scale_colour_gradient2(
      name = "FPR (%)",
      low = "dodgerblue3", mid = "white", high = "firebrick",
      midpoint = 50, limits = c(0, 100), na.value = "white"
    ) +
    scale_size_continuous(
      name = "FP count",
      range = c(3, 15),
      breaks = c(50, 500, 1000, 2000),
      labels = c("50", "0.5k", "1k", "2k")
    ) +
    facet_grid(cols = vars(width)) +
    labs(x = "SNP type", y = NULL, title = "S3A - False-positive heterozygous calls") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 10, height = 2.6, units = "in")
  message("Wrote ", outfile)
  invisible(df_summary)
}

plot_fnr_bubbles <- function(df1, outfile) {
  df_summary <- summarise_fnr(df1)
  p <- ggplot(df_summary, aes(x = new_SNP, y = 1)) +
    geom_point(aes(size = Hom_SNPs, colour = FNR)) +
    geom_text(aes(label = round(FNR, 1)), size = 3.5, colour = "black") +
    scale_colour_gradient2(
      name = "FNR (%)",
      low = "lightblue", mid = "white", high = "orange2",
      midpoint = 50, limits = c(0, 100), na.value = "white"
    ) +
    scale_size_continuous(
      name = "FN count",
      range = c(3, 15),
      breaks = seq(1e3, 1e4, 3e3),
      labels = paste0(seq(1, 10, 3), "k")
    ) +
    facet_grid(cols = vars(width)) +
    labs(x = "SNP type", y = NULL, title = "S3B - False-negative heterozygous calls") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 10, height = 2.6, units = "in")
  message("Wrote ", outfile)
  invisible(df_summary)
}

plot_baf_error_hist <- function(dt, outfile, title, which = c("het", "hom")) {
  which <- match.arg(which)
  d <- as.data.table(dt)
  d[, width := factor(
    as.character(width),
    levels = c("1", "2", "4"),
    labels = WIDTH_LABELS
  )]
  p <- ggplot(d, aes(x = BAF_diff_n)) +
    geom_histogram(
      fill = "white", colour = "grey40", binwidth = 0.02, linewidth = 0.2
    ) +
    facet_grid(cols = vars(width), scales = "free_y") +
    labs(
      x = "BAF_RRBS - BAF_WGS (normal)",
      y = "SNP count",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  ggsave(outfile, p, width = 7, height = 2.4, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

plot_fnr_fpr_tradeoff <- function(df_cov_thres2, outfile) {
  df_cov_thres <- rbindlist(df_cov_thres2, fill = TRUE)
  num_cols <- setdiff(names(df_cov_thres), c("thres", "min"))
  df_cov_thres <- df_cov_thres[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    .SDcols = num_cols,
    by = .(thres, min)
  ]
  df_cov_thres[, `:=`(
    thres = as.character(thres),
    min = as.character(min),
    FPR = FP / (FP + TN),
    FNR = FN / (FN + TP)
  )]
  p <- ggplot(df_cov_thres, aes(y = FNR, x = FPR, colour = thres)) +
    geom_point(aes(shape = min), size = 2) +
    geom_line(linewidth = 0.4) +
    scale_color_manual(
      name = "Heterozygosity\nboundaries",
      values = c(
        "0.1" = "grey10", "0.15" = "mediumpurple3", "0.2" = "cornflowerblue",
        "0.25" = "firebrick", "0.3" = "darkorange3"
      )
    ) +
    scale_shape_manual(
      name = "Minimum\ncoverage",
      values = c("10" = 16, "15" = 17, "20" = 4, "25" = 19, "30" = 8)
    ) +
    labs(
      x = "FPR", y = "FNR",
      title = "S3C - Heterozygous SNP calling performance"
    ) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  ggsave(outfile, p, width = 3.8, height = 2.8, units = "in")
  message("Wrote ", outfile)
  invisible(df_cov_thres)
}

plot_bias_sources <- function(df_plot, outfile) {
  d <- as.data.table(df_plot)
  setorder(d, type, flag1, flag2, flag4)
  d[, row_i := seq_len(.N)]
  p <- ggplot(d) +
    geom_tile(
      aes(x = 1, y = row_i),
      fill = ifelse(d$type == "Heterozygous", "grey50", "grey10")
    ) +
    geom_tile(
      aes(x = 2, y = row_i),
      fill = ifelse(d$flag1, "mediumpurple3", "grey99")
    ) +
    geom_tile(
      aes(x = 3, y = row_i),
      fill = ifelse(d$flag2, "orange3", "grey99")
    ) +
    geom_tile(
      aes(x = 4, y = row_i),
      fill = ifelse(!is.na(d$flag4) & d$flag4, "cornflowerblue", "grey99")
    ) +
    scale_x_continuous(
      name = "Source of bias",
      breaks = 1:4,
      labels = c(
        "RRBS\ngenotype",
        "At or near\nhet CCGG",
        "Short\nfragments",
        "Low ref\nread mq"
      )
    ) +
    labs(
      y = NULL,
      title = "S3D - TP / FN sources of bias (CRUK0062 normal)"
    ) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_text(size = 7),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 3.5, height = 5, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

################################################################################
# CLI
################################################################################

if (sys.nframe() == 0L) {

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figureS3",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory [default %default]"),
  make_option(c("--genotype_rdata"), type = "character",
              default = "~/rdeocd/txMeth/all_pilot_patients_genotype_TNR_FNR.RData",
              help = "Cached pilot genotype comparison object [default %default]"),
  make_option(c("--panels"), type = "character",
              default = "A,B,C,D",
              help = "Comma-separated panels [default %default]"),
  make_option(c("--thres"), type = "double",
              default = 0.1,
              help = "BAF heterozygosity boundary for A/B reclassify [default %default]"),
  make_option(c("--d_patient_index"), type = "integer",
              default = 2L,
              help = "df_plot2 index for S3D (1=CRUK0069,2=CRUK0062,3=CRUK0031) [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(work_dir)
output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

rdata <- path.expand(opt$genotype_rdata)
if (!file.exists(rdata)) {
  stop("Missing genotype RData: ", rdata)
}
message("Loading ", rdata)
e <- new.env(parent = emptyenv())
load(rdata, envir = e)
need <- c("df1_list", "df2_list", "df_cov_thres2", "df_plot2")
miss <- setdiff(need, ls(e))
if (length(miss) > 0) {
  stop("RData missing objects: ", paste(miss, collapse = ", "))
}

panels <- toupper(trimws(unlist(strsplit(opt$panels, ",", fixed = TRUE))))
thres <- as.numeric(opt$thres)

df1 <- reclassify_type(rbindlist(e$df1_list, fill = TRUE), thres = thres)
df2 <- reclassify_type(rbindlist(e$df2_list, fill = TRUE), thres = thres)

if ("A" %in% panels) {
  fpr_tab <- plot_fpr_bubbles(
    df2,
    file.path(output_dir, "S3A_FPR_by_SNP_type.pdf")
  )
  fwrite(fpr_tab, file.path(output_dir, "S3A_FPR_by_SNP_type.tsv"), sep = "\t")
  plot_baf_error_hist(
    df2,
    file.path(output_dir, "S3A_BAF_error_hist_WGS_hom.pdf"),
    title = "S3A - BAF error at WGS-homozygous SNPs",
    which = "hom"
  )
}

if ("B" %in% panels) {
  fnr_tab <- plot_fnr_bubbles(
    df1,
    file.path(output_dir, "S3B_FNR_by_SNP_type.pdf")
  )
  fwrite(fnr_tab, file.path(output_dir, "S3B_FNR_by_SNP_type.tsv"), sep = "\t")
  plot_baf_error_hist(
    df1,
    file.path(output_dir, "S3B_BAF_error_hist_WGS_het.pdf"),
    title = "S3B - BAF error at WGS-heterozygous SNPs",
    which = "het"
  )
}

if ("C" %in% panels) {
  cov_tab <- plot_fnr_fpr_tradeoff(
    e$df_cov_thres2,
    file.path(output_dir, "S3C_FNR_vs_FPR_thresholds.pdf")
  )
  fwrite(cov_tab, file.path(output_dir, "S3C_FNR_vs_FPR_thresholds.tsv"), sep = "\t")
}

if ("D" %in% panels) {
  idx <- as.integer(opt$d_patient_index)
  if (idx < 1L || idx > length(e$df_plot2)) {
    stop("--d_patient_index out of range for df_plot2")
  }
  # Notes final figure also rbindlist's all patients; legend cites CRUK0062
  plot_bias_sources(
    e$df_plot2[[idx]],
    file.path(output_dir, "S3D_TP_FN_sources_of_bias.pdf")
  )
}

message("Done.")

} # end main
