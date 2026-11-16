#' @description Correlate TRACERx RRBS–WES ASCN overlap with WES purity and
#' ploidy (WES as ground truth), and with RRBS–WES purity/ploidy disagreement.
#'
#' @author Lili Larose Cadieux
#' @return Prints Spearman / Pearson correlations and writes a summary table
#' and scatter plots under outputs/.

# Load library
x <- c("data.table", "ggplot2", "optparse")
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
# Source functions
################################################################################

#' Spearman and Pearson correlation with n and p-value
#' @param x Numeric vector
#' @param y Numeric vector
#' @return Named list with n, spearman_r, spearman_p, pearson_r, pearson_p
cor_stats <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  pe <- suppressWarnings(cor.test(x, y, method = "pearson"))
  list(
    n = length(x),
    spearman_r = unname(sp$estimate),
    spearman_p = sp$p.value,
    pearson_r = unname(pe$estimate),
    pearson_p = pe$p.value
  )
}

#' Format one correlation row for printing / tables
#' @param label Character label for the test
#' @param stats List from cor_stats()
#' @return data.table with one row
cor_row <- function(label, stats) {
  data.table(
    test = label,
    n = stats$n,
    spearman_r = round(stats$spearman_r, 3),
    spearman_p = stats$spearman_p,
    pearson_r = round(stats$pearson_r, 3),
    pearson_p = stats$pearson_p
  )
}

################################################################################
# parse command line arguments
################################################################################

option_list <- list(
  make_option(c("--wdir"),
              type = "character", default = NULL,
              help = "Working repo base directory (figure1H-K).")
)
opt <- parse_args(OptionParser(option_list = option_list))

# get working directory
work_dir <- ifelse(is.null(opt$wdir),
                   path.expand("~/data/CAMDAC_paper/figure1H-K"),
                   opt$wdir)
if (!file.exists(work_dir)) {
  stop(paste0("Working directory does not exist:", work_dir))
}
setwd(work_dir)

# set input and output directories
in_dir <- "inputs"
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# Load TRACERx RRBS–WES overlap table (WES purity / ploidy = ground truth)
################################################################################

load(file.path(in_dir, "dt_tracerx.RData"))
dt <- as.data.table(dt_tracerx)

dt[, perc_overlap := as.numeric(as.character(perc_copy_overlap))]
dt[, perc_overlap_1 := as.numeric(as.character(perc_copy_within_1))]

# disagreement vs WES ground truth
dt[, rho_diff := abs(rrbs_purity - wes_purity)]
dt[, psi_diff := abs(rrbs_ploidy - wes_ploidy)]

cat("N TRACERx RRBS–WES pairs:", nrow(dt), "\n\n")

################################################################################
# Correlations: perc_overlap vs WES purity / ploidy (ground truth)
################################################################################

res <- rbindlist(list(
  cor_row("perc_overlap ~ wes_purity",
          cor_stats(dt$perc_overlap, dt$wes_purity)),
  cor_row("perc_overlap ~ wes_ploidy",
          cor_stats(dt$perc_overlap, dt$wes_ploidy)),
  cor_row("perc_overlap_1 ~ wes_purity",
          cor_stats(dt$perc_overlap_1, dt$wes_purity)),
  cor_row("perc_overlap_1 ~ wes_ploidy",
          cor_stats(dt$perc_overlap_1, dt$wes_ploidy)),
  cor_row("rho_diff ~ wes_purity",
          cor_stats(dt$rho_diff, dt$wes_purity)),
  cor_row("rho_diff ~ wes_ploidy",
          cor_stats(dt$rho_diff, dt$wes_ploidy)),
  cor_row("psi_diff ~ wes_purity",
          cor_stats(dt$psi_diff, dt$wes_purity)),
  cor_row("psi_diff ~ wes_ploidy",
          cor_stats(dt$psi_diff, dt$wes_ploidy))
))

cat("=== Correlations (WES as ground truth) ===\n")
print(res)

out_tbl <- file.path(out_dir, "overlap_vs_wes_purity_ploidy_correlations.csv")
fwrite(res, out_tbl)
cat("\nWrote", out_tbl, "\n")

################################################################################
# Scatter plots: ASCN overlap vs WES purity / ploidy
################################################################################

theme_fig <- theme_classic() +
  theme(legend.position = "none")

p_rho <- ggplot(dt, aes(x = wes_purity, y = perc_overlap)) +
  geom_point(size = 2.5, alpha = 0.75, color = "cornflowerblue") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
  scale_x_continuous("WES purity (ground truth)",
                     breaks = seq(0.1, 0.9, 0.2), limits = c(0.1, 0.9)) +
  scale_y_continuous("RRBS-WES ASCN % overlap",
                     breaks = seq(0, 100, 20), limits = c(0, 100)) +
  theme_fig

ggsave(p_rho,
       filename = file.path(out_dir, "perc_overlap_vs_wes_purity.pdf"),
       device = "pdf", units = "in", width = 3.5, height = 3)

p_psi <- ggplot(dt, aes(x = wes_ploidy, y = perc_overlap)) +
  geom_point(size = 2.5, alpha = 0.75, color = "cornflowerblue") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
  scale_x_continuous("WES ploidy (ground truth)",
                     breaks = seq(2, 5, 1), limits = c(1.5, 5.5)) +
  scale_y_continuous("RRBS-WES ASCN % overlap",
                     breaks = seq(0, 100, 20), limits = c(0, 100)) +
  theme_fig

ggsave(p_psi,
       filename = file.path(out_dir, "perc_overlap_vs_wes_ploidy.pdf"),
       device = "pdf", units = "in", width = 3.5, height = 3)

# disagreement panels vs WES ground truth
p_rho_diff <- ggplot(dt, aes(x = wes_purity, y = rho_diff)) +
  geom_point(size = 2.5, alpha = 0.75, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
  scale_x_continuous("WES purity (ground truth)",
                     breaks = seq(0.1, 0.9, 0.2), limits = c(0.1, 0.9)) +
  scale_y_continuous("|RRBS - WES| purity") +
  theme_fig

ggsave(p_rho_diff,
       filename = file.path(out_dir, "rho_diff_vs_wes_purity.pdf"),
       device = "pdf", units = "in", width = 3.5, height = 3)

p_psi_diff <- ggplot(dt, aes(x = wes_ploidy, y = psi_diff)) +
  geom_point(size = 2.5, alpha = 0.75, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
  scale_x_continuous("WES ploidy (ground truth)",
                     breaks = seq(2, 5, 1), limits = c(1.5, 5.5)) +
  scale_y_continuous("|RRBS - WES| ploidy") +
  theme_fig

ggsave(p_psi_diff,
       filename = file.path(out_dir, "psi_diff_vs_wes_ploidy.pdf"),
       device = "pdf", units = "in", width = 3.5, height = 3)

# per-sample table for downstream use
out_samples <- file.path(out_dir, "rrbs_wes_overlap_vs_wes_gt.csv")
fwrite(dt[, .(
  patient_id, tx_id, rrbs_sample, wes_sample,
  wes_purity, wes_ploidy,
  rrbs_purity, rrbs_ploidy,
  rho_diff, psi_diff,
  perc_overlap, perc_overlap_1
)], out_samples)
cat("Wrote", out_samples, "\n")

# END
