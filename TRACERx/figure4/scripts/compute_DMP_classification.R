#!/usr/bin/env Rscript

#' Build per-sample ASM DMP classification tables for Figure 4 downstream plots.
#'
#' Writes \code{DMP_classification_{sample}.RData} (object \code{dt_DMP_classification})
#' under each \code{Allele_specific_analysis/{sample}/} folder.
#'
#' Requires prior phasing: \code{allele_specific_methylation_rates_{sample}.RData}
#' from \code{figure2P-S_S6F-I/scripts/Phase_methylation_at_SNVs.R}.
#'
#' Logic from \code{R_scripts/Notes_ASM_plots_16Jun2021.R} (\code{re_run=TRUE} block).

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
  optparse::make_option("--work_dir", type = "character",
                         default = "~/data/CAMDAC_paper/figure4"),
  optparse::make_option("--patients", type = "character", default = ""),
  optparse::make_option("--samples", type = "character", default = ""),
  optparse::make_option("--n_cores", type = "integer", default = 8L),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (!is.null(opt$patients) && nzchar(opt$patients)) {
  opt$patients <- paste(resolve_tx_patient_id(trimws(unlist(strsplit(opt$patients, ",", fixed = TRUE))), .tx_map), collapse = ",")
}

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

source(file.path(path.expand(opt$work_dir), "scripts", "asm_utils.R"))

split_arg <- function(x) {
  if (!nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

tx_dir <- path.expand(opt$tx_dir)
clinical <- file.path(
  tx_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
if (!file.exists(clinical)) stop("Clinical table not found: ", clinical)
patient_info <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
patient_info <- patient_info[!grepl("LN", Sample)]

patient_ids <- unique(patient_info$Patient)
patient_ids <- patient_ids[patient_ids != resolve_tx_patient_id("CRUK0047", .tx_map)]
req_pat <- split_arg(opt$patients)
req_smp <- split_arg(opt$samples)
if (length(req_pat)) patient_ids <- patient_ids[patient_ids %in% req_pat]

classify_one_sample <- function(patient_id, sample_id, rho) {
  path_patient <- file.path(tx_dir, patient_id)
  out_dir <- file.path(path_patient, "Allele_specific_analysis", sample_id)
  out_file <- file.path(out_dir, paste0("DMP_classification_", sample_id, ".RData"))
  if (file.exists(out_file) && !opt$overwrite) {
    return(list(status = "skip", file = out_file))
  }

  dmp_fls <- file.path(
    path_patient, "Methylation", sample_id, "re_segmented",
    "purified_tumour_DMP_calls.RData"
  )
  if (!file.exists(dmp_fls)) {
    dmp_fls <- file.path(
      path_patient, "Methylation", sample_id, "raw",
      "purified_tumour_DMP_calls.RData"
    )
  }
  if (!file.exists(dmp_fls)) {
    return(list(status = "missing_dmp", file = NA_character_))
  }

  phased_fls <- file.path(
    out_dir, paste0("allele_specific_methylation_rates_", sample_id, ".RData")
  )
  if (!file.exists(phased_fls)) {
    return(list(status = "missing_phased", file = NA_character_))
  }

  load(dmp_fls)
  df_mt <- data.table(df_mt)
  df_mt[, prob_DMP := ifelse(m_diff <= 0, 1 - prob_hypo, 1 - prob_hyper)]
  df_mt[, DMP_t := ifelse(
    m_diff < (-0.2) & prob_DMP > 0.99, "hypo",
    ifelse(m_diff > 0.2 & prob_DMP > 0.99, "hyper", NA)
  )]
  UM <- df_mt$m_n[df_mt$m_n_high < 0.2]
  M <- df_mt$m_n[df_mt$m_n_low > 0.7]
  df_mt <- df_mt[, .SD, .SDcols = c("chromosome", "start", "end", "prob_DMP", "DMP_t", "DMP_b")]
  setnames(df_mt, "chromosome", "chrom")
  setkeyv(df_mt, c("chrom", "start", "end"))

  load(phased_fls)
  df_phased_methylation[, DMP_t := NULL]
  dt <- merge(df_phased_methylation, df_mt, by = c("chrom", "start", "end"))
  rm(df_phased_methylation, df_mt)

  dt <- dt[!is.na(is_DM_t12)]
  dt <- dt[n1 != 0 & n2 != 0]
  dt[, multip := ifelse(n1 >= n2, paste(n1, n2, sep = "_"), paste(n2, n1, sep = "_"))]
  dt <- dt[(m_n_high < 0.15 | m_n_low > 0.8)]
  dt[, m_diff_tn := (m_t_corr - m_n)]

  alpha_n <- ceiling((dt$M_n + 1) / 2)
  beta_n <- ceiling((dt$UM_n + 1) / 2)
  alpha_b <- dt[, ceiling((M_1 + M_2 + 1) / 2)]
  beta_b <- dt[, ceiling((UM_1 + UM_2 + 1) / 2)]
  n <- nrow(dt)
  set.seed(1234)
  prob_hypo <- mcmapply(
    function(a_n, b_n, a_b, b_b) {
      out <- h(a_n, b_n, a_b, b_b)
      if (is.null(out)) NA else out
    },
    alpha_n, beta_n, alpha_b, beta_b,
    mc.cores = opt$n_cores
  )
  tmp_m_b <- alpha_b / (alpha_b + beta_b)
  tmp_m_n <- alpha_n / (alpha_n + beta_n)
  dt[, prob_DMP := ifelse(prob_hypo > 1, 1, ifelse(prob_hypo < 0, abs(prob_hypo), prob_hypo))]
  dt[, prob_DMP := ifelse(tmp_m_b - tmp_m_n > 0, 1 - prob_DMP, prob_DMP)]
  dt[, DMP_t := ifelse(
    m_diff_tn < (-0.2) & prob_DMP > 0.99, "hypo",
    ifelse(m_diff_tn > 0.2 & prob_DMP > 0.99, "hyper", NA)
  )]

  idx2 <- dt[, which(m_t_corr - m_n > 0.1)]
  selector <- dt[idx2, ifelse(m_t_corr - m_n >= 0, 0, 1)]
  selector[dt[idx2, m_t_corr - m_n == 0]] <- sample(
    c(0, 1), dt[idx2, sum(m_t_corr - m_n == 0)], replace = TRUE
  )
  prop_1 <- numeric(length(idx2))
  prop_1[selector == 0] <- M[sample(seq_along(M), sum(selector == 0), replace = TRUE)]
  prop_1[selector == 1] <- UM[sample(seq_along(UM), sum(selector == 1), replace = TRUE)]
  prop_2 <- numeric(length(idx2))
  prop_2[selector == 1] <- M[sample(seq_along(M), sum(selector == 1), replace = TRUE)]
  prop_2[selector == 0] <- UM[sample(seq_along(UM), sum(selector == 0), replace = TRUE)]

  dt[, cov_t := round2(
    ((total_counts_m_1 + total_counts_m_2) * rho * (n1 + n2)) /
      (rho * (n1 + n2) + (1 - rho) * (n1_n + n2_n)),
    digits = 0
  )]
  dt[, cov_n := total_counts_m_n_1 + total_counts_m_n_2]

  M_t <- unlist(mcmapply(function(cov, prop) rbinom(1, cov, prop),
                         cov = dt$cov_t[idx2], prop = prop_1, mc.cores = opt$n_cores))
  UM_t <- dt$cov_t[idx2] - M_t
  M_t <- M_t + 1
  UM_t <- UM_t + 1
  M_n2 <- unlist(mcmapply(function(cov, prop) rbinom(1, cov, prop),
                          cov = dt$cov_n[idx2], prop = prop_2, mc.cores = opt$n_cores))
  UM_n2 <- dt$cov_n[idx2] - M_n2
  M_n2 <- M_n2 + 1
  UM_n2 <- UM_n2 + 1
  dt[, ov_tn_sim := as.numeric(NA)]
  dt[idx2, ov_tn_sim := vapply(seq_along(idx2), function(i) {
    MCMC_overlap_beta(M_t[i], UM_t[i], M_n2[i], UM_n2[i])
  }, numeric(1))]
  dt[, pwr_tn := ov_tn_sim < 0.02]

  idx2 <- dt[, which(abs(m_t_1 - m_t_2) > 0.2)]
  selector <- dt[idx2, ifelse(m_t_1 - m_t_2 >= 0, 0, 1)]
  selector[dt[idx2, m_t_1 - m_t_2 == 0]] <- sample(
    c(0, 1), dt[idx2, sum(m_t_1 - m_t_2 == 0)], replace = TRUE
  )
  prop_1 <- numeric(length(idx2))
  prop_1[selector == 0] <- M[sample(seq_along(M), sum(selector == 0), replace = TRUE)]
  prop_1[selector == 1] <- UM[sample(seq_along(UM), sum(selector == 1), replace = TRUE)]
  prop_2 <- numeric(length(idx2))
  prop_2[selector == 1] <- M[sample(seq_along(M), sum(selector == 1), replace = TRUE)]
  prop_2[selector == 0] <- UM[sample(seq_along(UM), sum(selector == 0), replace = TRUE)]

  dt[, cov_t_1 := round2((total_counts_m_1 * rho * n1) / (rho * n1 + (1 - rho) * n1_n), digits = 0)]
  dt[, cov_t_2 := round2((total_counts_m_2 * rho * n2) / (rho * n2 + (1 - rho) * n2_n), digits = 0)]

  M_t1 <- unlist(mcmapply(function(cov, prop) rbinom(1, cov, prop),
                          cov = dt$cov_t_1[idx2], prop = prop_1, mc.cores = opt$n_cores))
  UM_t1 <- dt$cov_t_1[idx2] - M_t1
  M_t1 <- M_t1 + 1
  UM_t1 <- UM_t1 + 1
  M_t2 <- unlist(mcmapply(function(cov, prop) rbinom(1, cov, prop),
                          cov = dt$cov_t_2[idx2], prop = prop_2, mc.cores = opt$n_cores))
  UM_t2 <- dt$cov_t_2[idx2] - M_t2
  M_t2 <- M_t2 + 1
  UM_t2 <- UM_t2 + 1
  dt[, ov_t12_sim := as.numeric(NA)]
  dt[idx2, ov_t12_sim := vapply(seq_along(idx2), function(i) {
    MCMC_overlap_beta(M_t1[i], UM_t1[i], M_t2[i], UM_t2[i])
  }, numeric(1))]
  dt[, pwr_tt := ov_t12_sim < 0.02]

  dt[, is_DM_t12 := ((m_t_1_low > m_t_2_high) | (m_t_2_low > m_t_1_high)) & abs(m_diff_t12) > 0.4]
  dt[, is_IM_t := m_t_corr > 0.2 & m_t_corr < 0.7 & m_t_low > 0.15 & m_t_high < 0.75]

  dt[, flag := ifelse(
    !is.na(pwr_tt) & pwr_tt == TRUE & is_DM_t12 == TRUE & m_t_1 > m_t_2, "clonal-AS-DMP1",
    ifelse(
      !is.na(pwr_tt) & pwr_tt == TRUE & is_DM_t12 == TRUE & m_t_1 < m_t_2, "clonal-AS-DMP2",
      ifelse(
        !is.na(pwr_tn) & pwr_tn == TRUE & !is.na(DMP_t) & m_diff_t12 <= 0.4 & abs(m_t_corr - m_n) > 0.7,
        "clonal-bi-DMP",
        ifelse(
          !is.na(pwr_tn) & pwr_tn == TRUE & !is.na(DMP_t) & is_DM_t12 == FALSE & m_diff_t12 <= 0.4,
          "subclonal-DMP",
          ifelse(
            !is.na(pwr_tt) & !is.na(pwr_tn) & (pwr_tt == FALSE | pwr_tn == FALSE), "no-power",
            "non-DMP"
          )
        )
      )
    )
  )]
  dt[, flag := factor(
    flag,
    levels = c(
      "clonal-AS-DMP1", "clonal-AS-DMP2", "clonal-bi-DMP",
      "subclonal-DMP", "non-DMP", "no-power"
    ),
    ordered = TRUE
  )]

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dt_DMP_classification <- dt
  save(dt_DMP_classification, file = out_file)
  list(status = "ok", file = out_file, n = nrow(dt))
}

n_ok <- 0L
for (patient_id in patient_ids) {
  sample_ids <- patient_info[Patient == patient_id, unique(Sample)]
  if (length(req_smp)) sample_ids <- sample_ids[sample_ids %in% req_smp]
  for (sample_id in sample_ids) {
    rho <- as.numeric(patient_info[Patient == patient_id & Sample == sample_id, rho_reseg])
    cat(patient_id, sample_id, "... ")
    res <- tryCatch(
      classify_one_sample(patient_id, sample_id, rho),
      error = function(e) list(status = paste0("error:", conditionMessage(e)), file = NA)
    )
    cat(res$status, "\n")
    if (identical(res$status, "ok")) n_ok <- n_ok + 1L
  }
}
cat("Finished. Wrote or skipped", n_ok, "classification tables.\n")
