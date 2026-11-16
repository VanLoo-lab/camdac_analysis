#' @description Joint effect of WES purity and ploidy on TRACERx RRBS-WES
#' ASCN overlap: additive LM, partial/added-variable plots, and bivariate
#' color map with fitted iso-contours.
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
# parse command line arguments
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
# Load per-sample overlap table (WES = ground truth)
################################################################################

f_samples <- file.path(out_dir, "rrbs_wes_overlap_vs_wes_gt.csv")
if (!file.exists(f_samples)) {
  # fallback: rebuild from dt_tracerx
  load(file.path(in_dir, "dt_tracerx.RData"))
  dt <- as.data.table(dt_tracerx)
  dt[, perc_overlap := as.numeric(as.character(perc_copy_overlap))]
} else {
  dt <- fread(f_samples)
}

dt[, perc_overlap := as.numeric(perc_overlap)]
dt[, wes_purity := as.numeric(wes_purity)]
dt[, wes_ploidy := as.numeric(wes_ploidy)]

################################################################################
# Joint additive model
################################################################################

m_pur <- lm(perc_overlap ~ wes_purity, data = dt)
m_psi <- lm(perc_overlap ~ wes_ploidy, data = dt)
m_joint <- lm(perc_overlap ~ wes_purity + wes_ploidy, data = dt)
m_int <- lm(perc_overlap ~ wes_purity * wes_ploidy, data = dt)

dt[, fitted_joint := fitted(m_joint)]
dt[, residual := resid(m_joint)]

# z-scored predictors for comparable betas
dt[, z_rho := as.numeric(scale(wes_purity))]
dt[, z_psi := as.numeric(scale(wes_ploidy))]
m_z <- lm(perc_overlap ~ z_rho + z_psi, data = dt)

# partial / added-variable residuals
r_pur <- resid(lm(wes_purity ~ wes_ploidy, data = dt))
r_psi <- resid(lm(wes_ploidy ~ wes_purity, data = dt))
r_ov_pur <- resid(lm(perc_overlap ~ wes_ploidy, data = dt))
r_ov_psi <- resid(lm(perc_overlap ~ wes_purity, data = dt))
dt[, av_purity_x := r_pur]
dt[, av_purity_y := r_ov_pur]
dt[, av_ploidy_x := r_psi]
dt[, av_ploidy_y := r_ov_psi]

partial_r_pur <- cor(r_ov_pur, r_pur)
partial_p_pur <- cor.test(r_ov_pur, r_pur)$p.value
partial_r_psi <- cor(r_ov_psi, r_psi)
partial_p_psi <- cor.test(r_ov_psi, r_psi)$p.value

cat("N:", nrow(dt), "\n")
cat(sprintf("cor(wes_purity, wes_ploidy)=%.3f\n",
            cor(dt$wes_purity, dt$wes_ploidy)))
cat(sprintf("adj R2 purity=%.3f  ploidy=%.3f  joint=%.3f  +interaction=%.3f\n",
            summary(m_pur)$adj.r.squared,
            summary(m_psi)$adj.r.squared,
            summary(m_joint)$adj.r.squared,
            summary(m_int)$adj.r.squared))
cat("Joint coefficients (raw):\n")
print(summary(m_joint)$coefficients)
cat("Joint coefficients (z-scored):\n")
print(summary(m_z)$coefficients)
cat(sprintf("interaction p=%.3g\n",
            summary(m_int)$coefficients["wes_purity:wes_ploidy", "Pr(>|t|)"]))
cat(sprintf("partial r(overlap, purity | ploidy)=%.3f p=%.3g\n",
            partial_r_pur, partial_p_pur))
cat(sprintf("partial r(overlap, ploidy | purity)=%.3f p=%.3g\n",
            partial_r_psi, partial_p_psi))

################################################################################
# Summary table
################################################################################

sm <- summary(m_joint)$coefficients
res <- data.table(
  metric = c(
    "adj_R2_purity_only", "adj_R2_ploidy_only", "adj_R2_joint",
    "adj_R2_joint_interaction",
    "beta_wes_purity", "beta_wes_ploidy",
    "beta_z_rho", "beta_z_psi",
    "partial_r_purity_given_ploidy", "partial_r_ploidy_given_purity",
    "interaction_p", "cor_purity_ploidy"
  ),
  value = c(
    summary(m_pur)$adj.r.squared,
    summary(m_psi)$adj.r.squared,
    summary(m_joint)$adj.r.squared,
    summary(m_int)$adj.r.squared,
    sm["wes_purity", "Estimate"],
    sm["wes_ploidy", "Estimate"],
    summary(m_z)$coefficients["z_rho", "Estimate"],
    summary(m_z)$coefficients["z_psi", "Estimate"],
    partial_r_pur,
    partial_r_psi,
    summary(m_int)$coefficients["wes_purity:wes_ploidy", "Pr(>|t|)"],
    cor(dt$wes_purity, dt$wes_ploidy)
  ),
  p = c(
    summary(m_pur)$coefficients[2, 4],
    summary(m_psi)$coefficients[2, 4],
    NA_real_,
    summary(m_int)$coefficients["wes_purity:wes_ploidy", "Pr(>|t|)"],
    sm["wes_purity", "Pr(>|t|)"],
    sm["wes_ploidy", "Pr(>|t|)"],
    summary(m_z)$coefficients["z_rho", "Pr(>|t|)"],
    summary(m_z)$coefficients["z_psi", "Pr(>|t|)"],
    partial_p_pur,
    partial_p_psi,
    summary(m_int)$coefficients["wes_purity:wes_ploidy", "Pr(>|t|)"],
    cor.test(dt$wes_purity, dt$wes_ploidy)$p.value
  )
)
fwrite(res, file.path(out_dir, "perc_overlap_joint_model_summary.csv"))
fwrite(dt[, .(rrbs_sample, wes_purity, wes_ploidy, perc_overlap,
              fitted_joint, residual)],
       file.path(out_dir, "perc_overlap_joint_model_per_sample.csv"))

################################################################################
# Plots
################################################################################

theme_fig <- theme_classic(base_size = 11)

# 1) Recommended: bivariate purity x ploidy colored by overlap + fitted iso-lines
a <- coef(m_joint)[1]
b <- coef(m_joint)[2]
c_coef <- coef(m_joint)[3]
levels_k <- c(70, 80, 90, 95)
rho_seq <- seq(0.1, 0.9, length.out = 100)
cds <- rbindlist(lapply(levels_k, function(k) {
  data.table(
    wes_purity = rho_seq,
    wes_ploidy = (k - a - b * rho_seq) / c_coef,
    level = factor(k)
  )
}))
cds <- cds[wes_ploidy >= 1.5 & wes_ploidy <= 5.5]

p_bivar <- ggplot(dt, aes(x = wes_purity, y = wes_ploidy, color = perc_overlap)) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_line(data = cds, aes(x = wes_purity, y = wes_ploidy, group = level),
            inherit.aes = FALSE, color = "grey35", linetype = "dashed",
            linewidth = 0.4) +
  geom_text(
    data = cds[, .SD[which.min(abs(wes_purity - 0.72))], by = level],
    aes(x = wes_purity, y = wes_ploidy, label = paste0(level, "%")),
    inherit.aes = FALSE, size = 2.6, color = "grey25", vjust = -0.35
  ) +
  scale_color_gradientn(
    name = "ASCN %\noverlap",
    colours = c("red", "yellow", "blue"),
    values = scales::rescale(c(0, 50, 100)),
    limits = c(0, 100),
    breaks = c(0, 50, 100)
  ) +
  scale_x_continuous("WES purity (ground truth)",
                     breaks = seq(0.1, 0.9, 0.2), limits = c(0.1, 0.9)) +
  scale_y_continuous("WES ploidy (ground truth)",
                     breaks = seq(2, 5, 1), limits = c(1.5, 5.5)) +
  theme_fig

ggsave(p_bivar,
       filename = file.path(out_dir, "perc_overlap_bivariate_wes_purity_ploidy.pdf"),
       device = "pdf", units = "in", width = 4.6, height = 3.4)

# 2) Observed vs fitted joint model
p_fit <- ggplot(dt, aes(x = fitted_joint, y = perc_overlap)) +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(aes(color = wes_purity), size = 2.5, alpha = 0.85) +
  scale_color_gradient(name = "WES\npurity", low = "#fee5d9", high = "#a50f15") +
  scale_x_continuous(
    "Fitted ASCN % overlap (perc_overlap ~ wes_purity + wes_ploidy)",
    limits = c(50, 100)
  ) +
  scale_y_continuous("Observed ASCN % overlap", limits = c(20, 100)) +
  theme_fig

ggsave(p_fit,
       filename = file.path(out_dir, "perc_overlap_observed_vs_fitted_joint.pdf"),
       device = "pdf", units = "in", width = 4.6, height = 3.4)

# 3) Added-variable (partial residual) plots
p_av_pur <- ggplot(dt, aes(x = av_purity_x, y = av_purity_y)) +
  geom_point(size = 2.2, alpha = 0.75, color = "grey35") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick3", linewidth = 0.6) +
  xlab("WES purity | ploidy") +
  ylab("ASCN % overlap | ploidy") +
  theme_fig

p_av_psi <- ggplot(dt, aes(x = av_ploidy_x, y = av_ploidy_y)) +
  geom_point(size = 2.2, alpha = 0.75, color = "grey35") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick3", linewidth = 0.6) +
  xlab("WES ploidy | purity") +
  ylab("ASCN % overlap | purity") +
  theme_fig

pdf(file.path(out_dir, "perc_overlap_added_variable_wes.pdf"),
    width = 7, height = 3.2)
grid.arrange(p_av_pur, p_av_psi, ncol = 2)
dev.off()

cat("Wrote joint-model plots to", out_dir, "\n")

# END
