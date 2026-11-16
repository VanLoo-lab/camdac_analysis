#' @title Supplementary Figure S1A–F: mb / mn densities and HDI99
#' @author Lili Larose Cadieux
#' @description Normal and bulk-tumour methylation-fraction densities for
#' CRUK0062-R3 (CRUK0062-R3), HDI99 demonstration, and CN-stratified zooms for
#' CpGs confidently unmethylated in the matched normal.
#'
#' Panels:
#' * A — all CpGs: mn vs mb density
#' * B — simulated HDI99 at example unmethylated CpGs (beta intervals)
#' * C — confidently unmethylated in normal (HDI99 high <= 0.2)
#' * D — zoom of C (gained methylation across all copy numbers)
#' * E — as D, multiplicity 2+0
#' * F — as D, multiplicity 2+2
#'
#' Extracted from \code{Notes_plot_methylation_99_HDI.R}.

# Load libraries
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

#' Resolve a path that may be absolute, ~-expanded, or relative to base_dir
resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) {
    return(path)
  }
  file.path(base_dir, path)
}

#' Interval width for HDI optimisation (TeachingDemos / CAMDAC)
interval_width <- function(low_tail_pr, icdf_name, cred_mass, ...) {
  icdf_name(cred_mass + low_tail_pr, ...) - icdf_name(low_tail_pr, ...)
}

#' Highest density interval from an inverse CDF
#' @param icdf_name R inverse CDF (e.g. \code{qbeta})
#' @param cred_mass Desired HDI mass (default 0.99)
HDIofICDF <- function(icdf_name, cred_mass = 0.99, tol = 1e-4, ...) {
  incred_mass <- 1.0 - cred_mass
  opt_info <- optimize(
    f = interval_width,
    interval = c(0, incred_mass),
    icdf_name = icdf_name,
    cred_mass = cred_mass,
    tol = tol,
    ...
  )
  hdi_low_tail <- opt_info$minimum
  setNames(
    icdf_name(c(hdi_low_tail, cred_mass + hdi_low_tail), ...),
    nm = c("low", "high")
  )
}

#' Expected bulk mb bounds for clonal hypermethylation on n copies
#' (from unmethylated normal prior, OC/UC = 0.2 style bounds in Notes)
expected_hyper_bounds <- function(n_copies, purity, cn_total,
                                  low_m = 0.8, high_m = 1.0,
                                  normal_m_high = 0.2) {
  denom <- cn_total * purity + 2 * (1 - purity)
  list(
    low = (low_m * n_copies * purity) / denom,
    high = (high_m * n_copies * purity +
              normal_m_high * 2 * (1 - purity)) / denom
  )
}

#' Load CAMDAC purified tumour DMP table (\code{df_mt})
load_df_mt <- function(tx_dir, patient_id, sample_id, reseg = TRUE) {
  if (reseg) {
    path <- file.path(
      tx_dir, patient_id, "Methylation", sample_id, "re_segmented",
      "purified_tumour_DMP_calls.RData"
    )
  } else {
    path <- file.path(
      tx_dir, patient_id, "Methylation", sample_id,
      "purified_tumour_DMP_calls.RData"
    )
  }
  if (!file.exists(path)) {
    stop("Missing purified tumour DMP calls: ", path)
  }
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("df_mt", envir = e)) {
    stop("Object 'df_mt' not found in: ", path)
  }
  data.table(e$df_mt)
}

default_clinical_info <- function(tx_dir) {
  file.path(
    path.expand(tx_dir),
    paste0(
      "tx100clinicalinfo_RRBS_published_status_",
      "immune_clusters_purity_ploidy_CpG_cov.txt"
    )
  )
}

#' Resolve CLI patient id (CRUKID (or Patient)) -> list(ltx, cruk) for paths vs outputs
resolve_patient_ids <- function(patient_id, tx_dir, clinical_info = NULL) {
  if (is.null(clinical_info) || !nzchar(clinical_info)) {
    clinical_info <- default_clinical_info(tx_dir)
  }
  clinical_info <- path.expand(clinical_info)
  if (!file.exists(clinical_info)) {
    stop("Clinical table not found: ", clinical_info)
  }
  dt <- unique(fread(clinical_info)[, .(Patient, CRUKID)])
  if (grepl("^CRUK", patient_id)) {
    hit <- dt[CRUKID == patient_id]
    if (nrow(hit) == 0) stop("Unknown CRUKID: ", patient_id)
    return(list(ltx = hit$Patient[1], cruk = hit$CRUKID[1]))
  }
  hit <- dt[Patient == patient_id]
  if (nrow(hit) == 0) stop("Unknown Patient id: ", patient_id)
  list(ltx = hit$Patient[1], cruk = hit$CRUKID[1])
}

#' Try to read purity from TRACERx clinical table; NULL if unavailable
read_purity <- function(tx_dir, patient_id, sample_id,
                        clinical_info = NULL) {
  if (is.null(clinical_info) || !nzchar(clinical_info)) {
    clinical_info <- default_clinical_info(tx_dir)
  }
  clinical_info <- path.expand(clinical_info)
  if (!file.exists(clinical_info)) {
    return(NULL)
  }
  dt <- fread(clinical_info)
  ltx_id <- paste0(patient_id, "_", sample_id)
  if ("LTX_ID" %in% names(dt) && "rho_reseg" %in% names(dt)) {
    rho <- dt[LTX_ID == ltx_id, rho_reseg]
    if (length(rho) == 1 && is.finite(rho)) {
      return(as.numeric(rho))
    }
  }
  if (all(c("Patient", "Sample", "rho_reseg") %in% names(dt))) {
    rho <- dt[Patient == patient_id & Sample == sample_id, rho_reseg]
    if (length(rho) == 1 && is.finite(rho)) {
      return(as.numeric(rho))
    }
  }
  NULL
}

#' Shared density theme (match Notes)
theme_s1 <- function() {
  theme_classic()
}

#' Overlay expected clonal-hyper mb intervals for one ASCN state
add_expected_vlines <- function(p, purity, n_a, n_b) {
  cn <- n_a + n_b
  b_minor <- expected_hyper_bounds(n_b, purity, cn)
  b_major <- expected_hyper_bounds(n_a, purity, cn)
  b_all <- expected_hyper_bounds(cn, purity, cn)
  p +
    geom_vline(xintercept = b_minor$low, linetype = "dashed",
               colour = "orange3") +
    geom_vline(xintercept = b_minor$high, linetype = "dashed",
               colour = "orange3") +
    geom_vline(xintercept = b_major$low, linetype = "dashed",
               colour = "mediumpurple") +
    geom_vline(xintercept = b_major$high, linetype = "dashed",
               colour = "mediumpurple") +
    geom_vline(xintercept = b_all$low, linetype = "dashed",
               colour = "dodgerblue") +
    geom_vline(xintercept = b_all$high, linetype = "dashed",
               colour = "dodgerblue")
}

save_panel <- function(plot, outfile, width = 3.5, height = 3) {
  ggsave(
    filename = outfile, plot = plot, device = "pdf",
    width = width, height = height, units = "in"
  )
  message("Wrote ", outfile)
  invisible(plot)
}

################################################################################
# Panel plotters
################################################################################

#' S1B: HDI99 demonstration for example unmethylated beta posteriors
plot_s1b_hdi_demo <- function(outfile) {
  # Methylated counts in examples: Beta(n+1, 9n+1) ≈ mean 0.1 for n in vec
  # Plus sparse single-meth examples Beta(2, 29/69/99)
  vec <- c(1, 3, 7, 10)
  df <- data.table(low = 0, high = 0, point = 0, idx = 1:7)
  for (n in vec) {
    i <- which(vec == n)
    x <- HDIofICDF(
      icdf_name = qbeta, cred_mass = 0.99,
      shape1 = n + 1, shape2 = (n * 10 - n) + 1
    )
    df[i, `:=`(low = x[["low"]], high = x[["high"]], point = 0.1)]
    if (i != 1) {
      y <- HDIofICDF(
        icdf_name = qbeta, cred_mass = 0.99,
        shape1 = 1 + 1, shape2 = (n * 10 - 1) + 1
      )
      df[i + 3, `:=`(
        low = y[["low"]], high = y[["high"]],
        point = 1 / (n * 10)
      )]
    }
  }
  x_labs <- paste0(
    "Beta(", c(vec, rep(1, 3)), ",",
    c(9, 27, 63, 90, 29, 69, 99), ")"
  )
  p <- ggplot(df) +
    geom_point(aes(x = idx, y = point)) +
    geom_errorbar(aes(x = idx, ymin = low, ymax = high)) +
    scale_x_reverse(
      name = "", breaks = 1:7, labels = x_labs
    ) +
    geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey75") +
    ylim(c(0, 1)) +
    coord_flip() +
    theme_s1()
  save_panel(p, outfile)
}

#' Density of mn (dashed) and mb (solid grey)
plot_mb_mn_density <- function(dt, outfile, adjust = 1,
                               ylim_max = NULL, title = NULL) {
  p <- ggplot(dt, aes(y = after_stat(density))) +
    theme_s1() +
    geom_density(aes(x = m_b), adjust = adjust, colour = "grey60") +
    geom_density(aes(x = m_n), adjust = adjust, linetype = "dashed") +
    xlim(c(0, 1))
  if (!is.null(ylim_max)) {
    p <- p + coord_cartesian(ylim = c(0, ylim_max))
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  save_panel(p, outfile)
}

#' CN-stratified zoom with expected clonal-hyper vlines
plot_mb_mn_zoom_multi <- function(dt, outfile, multi_keep, purity,
                                  n_a, n_b, adjust = 5, ylim_max = 0.7) {
  if (length(multi_keep) == 1) {
    dt_sub <- dt[multi == multi_keep]
  } else {
    dt_sub <- dt[multi %in% multi_keep]
  }
  p <- ggplot(dt_sub, aes(y = after_stat(density))) +
    theme_s1() +
    geom_density(aes(x = m_b), adjust = adjust, colour = "grey60") +
    geom_density(aes(x = m_n), adjust = adjust, linetype = "dashed") +
    xlim(c(0, 1)) +
    coord_cartesian(ylim = c(0, ylim_max))
  p <- add_expected_vlines(p, purity = purity, n_a = n_a, n_b = n_b)
  save_panel(p, outfile)
}

################################################################################
# Parse command line arguments
################################################################################

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figureS1A-F",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory relative to --work_dir or absolute [default %default]"),
  make_option(c("--tx_dir"), type = "character",
              default = "~/rdeocd/txMeth",
              help = "TRACERx methylation parent directory [default %default]"),
  make_option(c("--patient_id"), type = "character",
              default = "CRUK0062",
              help = "Patient ID as CRUKID (or Patient) (outputs use CRUKID) [default %default]"),
  make_option(c("--sample_id"), type = "character",
              default = "R3",
              help = "Tumour region / sample ID [default %default]"),
  make_option(c("--reseg"), type = "logical",
              default = TRUE,
              help = "Use re-segmented CAMDAC outputs [default %default]"),
  make_option(c("--purity"), type = "double",
              default = NA,
              help = paste0(
                "Tumour purity for expected-peak vlines (E–F). ",
                "Default: rho_reseg from clinical table, else 0.73 (Notes)"
              )),
  make_option(c("--clinical_info"), type = "character",
              default = NULL,
              help = "Optional clinical TSV override for purity lookup"),
  make_option(c("--panels"), type = "character",
              default = "A,B,C,D,E,F",
              help = "Comma-separated panels to write [default %default]"),
  make_option(c("--unmeth_hdi_max"), type = "double",
              default = 0.2,
              help = "m_n_high threshold for confidently unmethylated [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!is.null(opt$patient_id) && nzchar(opt$patient_id)) {
  opt$patient_id <- paste(resolve_tx_patient_id(trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE))), .tx_map), collapse = ",")
}

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}
setwd(work_dir)

output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tx_dir <- path.expand(opt$tx_dir)
ids <- resolve_patient_ids(opt$patient_id, tx_dir, opt$clinical_info)
patient_id <- ids$ltx
cruk_id <- ids$cruk
sample_id <- opt$sample_id
reseg <- isTRUE(opt$reseg)
panels <- toupper(trimws(unlist(strsplit(opt$panels, ",", fixed = TRUE))))
unmeth_max <- opt$unmeth_hdi_max
message("Patient ", patient_id, " -> output id ", cruk_id)

purity <- opt$purity
if (is.na(purity)) {
  purity <- read_purity(
    tx_dir, patient_id, sample_id,
    clinical_info = opt$clinical_info
  )
  if (is.null(purity)) {
    purity <- 0.73
    message("Using Notes default purity = ", purity)
  } else {
    message("Purity from clinical table: ", purity)
  }
}

prefix <- paste0(cruk_id, "_", sample_id)

################################################################################
# Generate panels
################################################################################

# Panel B needs no methylation data
if ("B" %in% panels) {
  plot_s1b_hdi_demo(
    file.path(output_dir, paste0(prefix, "_S1B_HDI_pointrange.pdf"))
  )
}

need_mt <- any(panels %in% c("A", "C", "D", "E", "F"))
if (need_mt) {
  message("Loading df_mt: ", patient_id, " ", sample_id)
  df_mt <- load_df_mt(tx_dir, patient_id, sample_id, reseg = reseg)
  req <- c("m_b", "m_n", "m_n_high")
  missing <- setdiff(req, names(df_mt))
  if (length(missing) > 0) {
    stop("df_mt missing columns: ", paste(missing, collapse = ", "))
  }
  if (!"multi" %in% names(df_mt) && any(panels %in% c("E", "F"))) {
    stop("df_mt missing column 'multi' required for panels E–F")
  }
  df_um <- df_mt[m_n_high <= unmeth_max]
}

if ("A" %in% panels) {
  # All CpGs — normal vs bulk density
  plot_mb_mn_density(
    df_mt,
    file.path(output_dir, paste0(prefix, "_S1A_all_CpGs_mb_mn.pdf")),
    adjust = 1
  )
}

if ("C" %in% panels) {
  plot_mb_mn_density(
    df_um,
    file.path(output_dir, paste0(prefix, "_S1C_normal_unmeth_mb_mn.pdf")),
    adjust = 5
  )
}

if ("D" %in% panels) {
  # Zoom highlighting gained methylation (all copy numbers)
  plot_mb_mn_density(
    df_um,
    file.path(output_dir, paste0(prefix, "_S1D_normal_unmeth_mb_mn_zoom.pdf")),
    adjust = 5, ylim_max = 0.75
  )
}

if ("E" %in% panels) {
  # 2+0 — fix Notes bug: compute nA/nB intervals before plotting
  plot_mb_mn_zoom_multi(
    df_um,
    file.path(output_dir, paste0(prefix, "_S1E_normal_unmeth_mb_mn_zoom_2plus0.pdf")),
    multi_keep = "2+0", purity = purity, n_a = 2, n_b = 0,
    adjust = 5, ylim_max = 0.7
  )
}

if ("F" %in% panels) {
  plot_mb_mn_zoom_multi(
    df_um,
    file.path(output_dir, paste0(prefix, "_S1F_normal_unmeth_mb_mn_zoom_2plus2.pdf")),
    multi_keep = "2+2", purity = purity, n_a = 2, n_b = 2,
    adjust = 5, ylim_max = 0.7
  )
}

message("Done.")
