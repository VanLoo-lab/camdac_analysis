#' @description Joint effect of DNA-seq purity and ploidy on BS-seq vs DNA-seq
#' ASCN overlap across matched-normal cohorts (TRACERx RRBS-WES, CRC RRBS-WGS,
#' CRC WGBS-WGS) and, separately, the mCRPC WGBS-WGS panel-of-normals set.
#'
#' @author Lili Larose Cadieux
#' @return Writes summary tables and PDFs under outputs/.

# Load library
x <- c("data.table", "ggplot2", "gridExtra", "optparse")
invisible(lapply(x, function(y)
  suppressWarnings(
    suppressMessages(
      library(y,
              character.only = TRUE,
              quietly = TRUE,
              warn.conflicts = FALSE))
    )
  )
)
rm(x)

################################################################################
# Helpers
################################################################################

fit_joint <- function(dt, label) {
  stopifnot(all(c("perc_overlap", "dnaseq_rho", "dnaseq_psi") %in% names(dt)))
  dt <- dt[complete.cases(perc_overlap, dnaseq_rho, dnaseq_psi)]
  m_pur <- lm(perc_overlap ~ dnaseq_rho, data = dt)
  m_psi <- lm(perc_overlap ~ dnaseq_psi, data = dt)
  m_joint <- lm(perc_overlap ~ dnaseq_rho + dnaseq_psi, data = dt)
  m_int <- lm(perc_overlap ~ dnaseq_rho * dnaseq_psi, data = dt)

  dt[, z_rho := as.numeric(scale(dnaseq_rho))]
  dt[, z_psi := as.numeric(scale(dnaseq_psi))]
  m_z <- lm(perc_overlap ~ z_rho + z_psi, data = dt)

  r_pur <- resid(lm(dnaseq_rho ~ dnaseq_psi, data = dt))
  r_psi <- resid(lm(dnaseq_psi ~ dnaseq_rho, data = dt))
  r_ov_pur <- resid(lm(perc_overlap ~ dnaseq_psi, data = dt))
  r_ov_psi <- resid(lm(perc_overlap ~ dnaseq_rho, data = dt))

  sm <- summary(m_joint)$coefficients
  smz <- summary(m_z)$coefficients
  smi <- summary(m_int)$coefficients

  summary_dt <- data.table(
    cohort = label,
    n = nrow(dt),
    metric = c(
      "adj_R2_purity_only", "adj_R2_ploidy_only", "adj_R2_joint",
      "adj_R2_joint_interaction",
      "beta_dnaseq_rho", "beta_dnaseq_psi",
      "beta_z_rho", "beta_z_psi",
      "partial_r_purity_given_ploidy", "partial_r_ploidy_given_purity",
      "spearman_overlap_purity", "spearman_overlap_ploidy",
      "interaction_p", "cor_purity_ploidy"
    ),
    value = c(
      summary(m_pur)$adj.r.squared,
      summary(m_psi)$adj.r.squared,
      summary(m_joint)$adj.r.squared,
      summary(m_int)$adj.r.squared,
      sm["dnaseq_rho", "Estimate"],
      sm["dnaseq_psi", "Estimate"],
      smz["z_rho", "Estimate"],
      smz["z_psi", "Estimate"],
      cor(r_ov_pur, r_pur),
      cor(r_ov_psi, r_psi),
      cor(dt$perc_overlap, dt$dnaseq_rho, method = "spearman"),
      cor(dt$perc_overlap, dt$dnaseq_psi, method = "spearman"),
      smi["dnaseq_rho:dnaseq_psi", "Pr(>|t|)"],
      cor(dt$dnaseq_rho, dt$dnaseq_psi)
    ),
    p = c(
      summary(m_pur)$coefficients[2, 4],
      summary(m_psi)$coefficients[2, 4],
      NA_real_,
      smi["dnaseq_rho:dnaseq_psi", "Pr(>|t|)"],
      sm["dnaseq_rho", "Pr(>|t|)"],
      sm["dnaseq_psi", "Pr(>|t|)"],
      smz["z_rho", "Pr(>|t|)"],
      smz["z_psi", "Pr(>|t|)"],
      cor.test(r_ov_pur, r_pur)$p.value,
      cor.test(r_ov_psi, r_psi)$p.value,
      suppressWarnings(cor.test(dt$perc_overlap, dt$dnaseq_rho,
                                method = "spearman")$p.value),
      suppressWarnings(cor.test(dt$perc_overlap, dt$dnaseq_psi,
                                method = "spearman")$p.value),
      smi["dnaseq_rho:dnaseq_psi", "Pr(>|t|)"],
      cor.test(dt$dnaseq_rho, dt$dnaseq_psi)$p.value
    )
  )

  dt[, fitted_joint := fitted(m_joint)]
  list(dt = dt, model = m_joint, summary = summary_dt)
}

plot_bivariate <- function(dt, m_joint, out_file, title_lab = NULL) {
  a <- coef(m_joint)[1]
  b <- coef(m_joint)[2]
  c_coef <- coef(m_joint)[3]
  levels_k <- c(70, 80, 90, 95)
  rho_seq <- seq(0.05, 0.95, length.out = 100)
  cds <- rbindlist(lapply(levels_k, function(k) {
    data.table(
      dnaseq_rho = rho_seq,
      dnaseq_psi = (k - a - b * rho_seq) / c_coef,
      level = factor(k)
    )
  }))
  cds <- cds[dnaseq_psi >= 1.2 & dnaseq_psi <= 5.8]

  p <- ggplot(dt, aes(x = dnaseq_rho, y = dnaseq_psi, color = perc_overlap)) +
    geom_point(aes(shape = comparison), size = 2.6, alpha = 0.9) +
    geom_line(data = cds, aes(x = dnaseq_rho, y = dnaseq_psi, group = level),
              inherit.aes = FALSE, color = "grey35", linetype = "dashed",
              linewidth = 0.4) +
    scale_color_gradientn(
      name = "ASCN %\noverlap",
      colours = c("red", "yellow", "blue"),
      values = scales::rescale(c(0, 50, 100)),
      limits = c(0, 100),
      breaks = c(0, 50, 100)
    ) +
    scale_shape_manual(
      name = "comparison",
      values = c(
        "RRBS vs WES" = 15,
        "RRBS vs WGS" = 16,
        "WGBS vs WGS" = 18
      )
    ) +
    scale_x_continuous("DNA-seq purity (ground truth)",
                       breaks = seq(0.1, 0.9, 0.2)) +
    scale_y_continuous("DNA-seq ploidy (ground truth)",
                       breaks = seq(1, 6, 1)) +
    theme_classic(base_size = 11)

  if (!is.null(title_lab)) {
    p <- p + ggtitle(title_lab) +
      theme(plot.title = element_text(size = 10, face = "plain"))
  }

  ggsave(p, filename = out_file, device = "pdf",
         units = "in", width = 5.2, height = 3.6)
  invisible(p)
}

################################################################################
# parse args / paths
################################################################################

option_list <- list(
  make_option(c("--wdir"),
              type = "character", default = NULL,
              help = "Working repo base directory (figure1H-K).")
)
opt <- parse_args(OptionParser(option_list = option_list))

work_dir <- ifelse(is.null(opt$wdir),
                   path.expand("~/data/CAMDAC_paper/figure1H-K"),
                   opt$wdir)
if (!file.exists(work_dir)) {
  stop(paste0("Working directory does not exist:", work_dir))
}
setwd(work_dir)

in_dir <- "inputs"
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# Build matched-normal table (Fig: bsseq_vs_dnaseq_ascn_overlap.pdf)
################################################################################

load(file.path(in_dir, "dt_tracerx.RData"))
dt_tracerx <- as.data.table(dt_tracerx)

dt_matched_list <- list(
  dt_tracerx[, .(
    sample_id = rrbs_sample,
    dnaseq_rho = as.numeric(wes_purity),
    dnaseq_psi = as.numeric(wes_ploidy),
    bsseq_rho = as.numeric(rrbs_purity),
    bsseq_psi = as.numeric(rrbs_ploidy),
    perc_overlap = as.numeric(as.character(perc_copy_overlap)),
    perc_overlap_1 = as.numeric(as.character(perc_copy_within_1)),
    comparison = "RRBS vs WES",
    normal_type = "matched",
    cohort = "TRACERx"
  )],
  fread(file.path(in_dir, "CRC_RRBS_WGS_overlap.txt"))[, .(
    sample_id,
    dnaseq_rho = as.numeric(purity_a),
    dnaseq_psi = as.numeric(ploidy_a),
    bsseq_rho = as.numeric(purity_b),
    bsseq_psi = as.numeric(ploidy_b),
    perc_overlap = as.numeric(perc_overlap),
    perc_overlap_1 = as.numeric(perc_overlap_1),
    comparison = "RRBS vs WGS",
    normal_type = "matched",
    cohort = "CRC"
  )],
  fread(file.path(in_dir, "CRC_WGBS_WGS_overlap.txt"))[, .(
    sample_id,
    dnaseq_rho = as.numeric(purity_a),
    dnaseq_psi = as.numeric(ploidy_a),
    bsseq_rho = as.numeric(purity_b),
    bsseq_psi = as.numeric(ploidy_b),
    perc_overlap = as.numeric(perc_overlap),
    perc_overlap_1 = as.numeric(perc_overlap_1),
    comparison = "WGBS vs WGS",
    normal_type = "matched",
    cohort = "CRC"
  )]
)

# optional TRACERx RRBS-WGS overlaps (from rrbs_wgs_intersections.R)
tx_wgs_f <- file.path(in_dir, "TRACERx_RRBS_WGS_overlap.txt")
if (file.exists(tx_wgs_f)) {
  dt_matched_list[[length(dt_matched_list) + 1]] <-
    fread(tx_wgs_f)[, .(
      sample_id,
      dnaseq_rho = as.numeric(purity_a),
      dnaseq_psi = as.numeric(ploidy_a),
      bsseq_rho = as.numeric(purity_b),
      bsseq_psi = as.numeric(ploidy_b),
      perc_overlap = as.numeric(perc_overlap),
      perc_overlap_1 = as.numeric(perc_overlap_1),
      comparison = "RRBS vs WGS",
      normal_type = "matched",
      cohort = "TRACERx"
    )]
}

dt_matched <- rbindlist(dt_matched_list, use.names = TRUE)

################################################################################
# Build panel-of-normals table (Fig: ..._panel_of_normals.pdf)
################################################################################

dt_pon <- fread(file.path(in_dir, "mcrpc_wgbs_wgs_overlap.csv"),
                sep = ",",
                colClasses = c(perc_overlap = "numeric"))
# same exclusion as figure1H-K.R
dt_pon <- dt_pon[sample_id != "DTB-255"]
dt_pon <- dt_pon[, .(
  sample_id,
  dnaseq_rho = as.numeric(purity_a),
  dnaseq_psi = as.numeric(ploidy_a),
  bsseq_rho = as.numeric(purity_b),
  bsseq_psi = as.numeric(ploidy_b),
  perc_overlap = as.numeric(perc_overlap),
  perc_overlap_1 = as.numeric(perc_overlap_1),
  comparison = "WGBS vs WGS",
  normal_type = "panel_of_normals"
)]

cat("Matched-normal n by comparison / cohort:\n")
print(dt_matched[, .N, by = .(comparison, cohort)])
cat("Panel-of-normals n:", nrow(dt_pon), "\n\n")

################################################################################
# Fit models
################################################################################

# overall matched + per-comparison matched + panel
fit_all_matched <- fit_joint(copy(dt_matched), "matched_all")
fit_tx_wes <- fit_joint(dt_matched[comparison == "RRBS vs WES"],
                        "matched_RRBS_vs_WES")
fit_rrbs_wgs <- fit_joint(dt_matched[comparison == "RRBS vs WGS"],
                          "matched_RRBS_vs_WGS")
fit_tx_wgs <- fit_joint(
  dt_matched[comparison == "RRBS vs WGS" & cohort == "TRACERx"],
  "matched_TRACERx_RRBS_vs_WGS"
)
fit_crc_wgbs <- fit_joint(dt_matched[comparison == "WGBS vs WGS"],
                          "matched_WGBS_vs_WGS")
fit_pon <- fit_joint(copy(dt_pon), "panel_of_normals_WGBS_vs_WGS")

res <- rbindlist(list(
  fit_all_matched$summary,
  fit_tx_wes$summary,
  fit_rrbs_wgs$summary,
  fit_tx_wgs$summary,
  fit_crc_wgbs$summary,
  fit_pon$summary
))

# compact key comparison table
key_metrics <- c(
  "adj_R2_joint", "adj_R2_purity_only", "adj_R2_ploidy_only",
  "beta_dnaseq_rho", "beta_dnaseq_psi",
  "spearman_overlap_purity", "spearman_overlap_ploidy",
  "partial_r_purity_given_ploidy", "partial_r_ploidy_given_purity"
)
res_key <- res[metric %in% key_metrics]
res_wide <- dcast(res_key, cohort + n ~ metric, value.var = "value")

cat("=== Key model summary (adj R2 / betas) ===\n")
print(res_wide[, .(
  cohort, n,
  adj_R2_joint = round(adj_R2_joint, 3),
  adj_R2_purity = round(adj_R2_purity_only, 3),
  adj_R2_ploidy = round(adj_R2_ploidy_only, 3),
  beta_rho = round(beta_dnaseq_rho, 2),
  beta_psi = round(beta_dnaseq_psi, 2),
  spearman_rho = round(spearman_overlap_purity, 3),
  spearman_psi = round(spearman_overlap_ploidy, 3)
)])

fwrite(res, file.path(out_dir, "perc_overlap_joint_all_cohorts_summary.csv"))
fwrite(res_wide, file.path(out_dir, "perc_overlap_joint_all_cohorts_key.csv"))

fwrite(fit_all_matched$dt[, .(
  sample_id, comparison, normal_type,
  dnaseq_rho, dnaseq_psi, bsseq_rho, bsseq_psi,
  perc_overlap, perc_overlap_1, fitted_joint
)], file.path(out_dir, "perc_overlap_joint_matched_per_sample.csv"))

fwrite(fit_pon$dt[, .(
  sample_id, comparison, normal_type,
  dnaseq_rho, dnaseq_psi, bsseq_rho, bsseq_psi,
  perc_overlap, perc_overlap_1, fitted_joint
)], file.path(out_dir, "perc_overlap_joint_panel_of_normals_per_sample.csv"))

################################################################################
# Plots
################################################################################

r2_m <- round(summary(fit_all_matched$model)$adj.r.squared, 3)
r2_p <- round(summary(fit_pon$model)$adj.r.squared, 3)

plot_bivariate(
  fit_all_matched$dt,
  fit_all_matched$model,
  file.path(out_dir, "perc_overlap_bivariate_matched_normals.pdf"),
    title_lab = sprintf(
    "Matched normals (TRACERx RRBS-WES/WGS + CRC RRBS/WGBS-WGS); joint adj R2=%.3f",
    r2_m
  )
)

plot_bivariate(
  fit_pon$dt,
  fit_pon$model,
  file.path(out_dir, "perc_overlap_bivariate_panel_of_normals.pdf"),
  title_lab = sprintf(
    "Panel of normals (mCRPC WGBS-WGS); joint adj R2=%.3f",
    r2_p
  )
)

# side-by-side observed vs fitted for matched vs PoN
dt_fit_cmp <- rbind(
  fit_all_matched$dt[, .(perc_overlap, fitted_joint, set = "Matched normals")],
  fit_pon$dt[, .(perc_overlap, fitted_joint, set = "Panel of normals")]
)
p_fit <- ggplot(dt_fit_cmp, aes(x = fitted_joint, y = perc_overlap)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(size = 2.2, alpha = 0.75, color = "grey30") +
  facet_wrap(~ set) +
  scale_x_continuous("Fitted ASCN % overlap\n(perc_overlap ~ DNA-seq purity + ploidy)") +
  scale_y_continuous("Observed ASCN % overlap", limits = c(0, 100)) +
  theme_classic(base_size = 11)

ggsave(p_fit,
       filename = file.path(out_dir, "perc_overlap_observed_vs_fitted_matched_vs_pon.pdf"),
       device = "pdf", units = "in", width = 7, height = 3.4)

# bar comparison of adj R2
r2_bars <- res[metric == "adj_R2_joint", .(
  cohort, n, adj_R2 = value
)]
r2_bars[, label := paste0(cohort, " (n=", n, ")")]
p_r2 <- ggplot(r2_bars, aes(x = reorder(label, adj_R2), y = adj_R2, fill = cohort)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous("Joint model adj. R2\n(perc_overlap ~ DNA-seq purity + ploidy)") +
  xlab(NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")

ggsave(p_r2,
       filename = file.path(out_dir, "perc_overlap_joint_adjR2_by_cohort.pdf"),
       device = "pdf", units = "in", width = 6.5, height = 3.2)

cat("\nWrote matched-normal and panel-of-normals joint-model outputs to",
    out_dir, "\n")

# END
