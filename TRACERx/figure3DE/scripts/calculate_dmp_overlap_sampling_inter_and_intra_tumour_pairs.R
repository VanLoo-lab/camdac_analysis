#' @description Sample CpG-wise tumour–tumour pairs (pure tumour mt and bulk mb
#' intervals) from genome-wide CAMDAC results, classify purity groups by
#' \code{rho_reseg}, and compute raw / corrected DMP-style overlap calls for
#' \code{inter}-patient versus \code{intra}-patient designs (v2 logic from
#' \code{Notes_mt_mt_vs_mb_mb_DMPs_figure5FG.R}). Produces Euler Venn PDFs per
#' purity comparison label.
#'
#' @details Purity groups: \code{group_rho == 1} (low, rho \< 0.3) and
#' \code{group_rho == 4} (high, rho \> 0.58). Copy-number stratification is
#' implemented in the Notes script but the CN loops are left commented there;
#' this script follows the same default with \code{CN1 <- NULL} (purity-only
#' stratification). Output locations match the Notes block: Euler PDFs under
#' \verb{~/rdeocd/txMeth_analysis/mb_mb_vs_mt_mt_DMPs/venn_plots/} and
#' \verb{power_real_results_list.RData} under \verb{~/rdeocd/txMeth/}. The
#' workspace cleanup \code{rm(list = ls(pattern = "CN|gr"))} from the Notes is
#' omitted here because it is unsafe and does not change saved results.
#'
#' @author Lili Larose Cadieux

################################################################################
# Load libraries
################################################################################

x <- c("data.table", "fst", "DescTools", "eulerr")

invisible(lapply(x, function(y) {
  suppressWarnings(
    suppressMessages(
      library(y, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
  )
}))
rm(x)

################################################################################
# Paths
################################################################################

parent_dir <- path.expand("~/rdeocd/txMeth")
fst_path <- file.path(parent_dir, "CpG_wise_CAMDAC_results_noNAs.fst")
clinical_path <- file.path(
  parent_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)

out_dir <- path.expand("~/data/CAMDAC_paper/figure3DE/outputs/dmp_overlap_sampling")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rdata_out <- file.path(out_dir, "power_real_results_list.RData")

################################################################################
# Clinical metadata and purity groups
################################################################################

patient_info <- as.data.table(
  read.delim(clinical_path, header = TRUE, sep = "\t")
)
patient_info[, group_rho := ifelse(
  rho_reseg < 0.3, 1,
  ifelse(rho_reseg > 0.58, 4, 0)
)]
table(patient_info$group_rho)
patient_ids <- unique(patient_info$Patient)

################################################################################
# CpG-wise CAMDAC results (autosomes only; drop chr X)
################################################################################

dt <- read.fst(
  fst_path,
  columns = c(
    "ID", "patient", "sample", "chromosome", "start", "end",
    "m_t_low", "m_t_high", "m_t",
    "m_b_low", "m_b_high", "m_b",
    "CN"
  ),
  as.data.table = TRUE
)
dt <- dt[chromosome != "X"]
dt[, LTXID := paste(patient, sample, sep = "_")]
dt[, group_rho := patient_info$group_rho[match(LTXID, patient_info$LTX_ID)]]
dt[, rho := patient_info$rho_reseg[match(LTXID, patient_info$LTX_ID)]]
dt <- dt[, .SD, .SDcols = c(
  "ID", "LTXID", "patient", "chromosome", "start", "end",
  "m_t_low", "m_t_high", "m_t", "m_b_low", "m_b_high", "m_b",
  "CN", "group_rho", "rho"
)]

# Restrict to low- or high-purity samples only
ids <- patient_info$LTX_ID[patient_info$group_rho == 1L | patient_info$group_rho == 4L]
dt <- dt[LTXID %in% ids]
rm(ids)

################################################################################
# Static parameters
################################################################################

set.seed(1234)
flags <- c("intra", "inter")
CNs <- seq(1, 6, 1)
CN1 <- NULL
grs <- c(1, 4)
reseg <- TRUE
effect_size <- 0.2

################################################################################
# Inter- vs intra-tumour sampling and DMP overlap
################################################################################

for (flag in flags) {
  dt_list <- list()

  for (gr1 in grs) {
    for (gr2 in grs) {
      if (!is.null(CN1)) {
        print(paste0(
          "copy number ", CN1, " vs. ", CN2,
          " and purity ", gr1, " vs. ", gr2
        ))
      } else {
        print(paste0("Purity groups ", gr1, " vs. ", gr2))
      }

      ids <- patient_info$LTX_ID[
        patient_info$group_rho == gr1 | patient_info$group_rho == gr2
      ]
      dt_tmp <- dt[LTXID %in% ids]
      rm(ids)

      if (!is.null(CN1)) {
        if (CN2 != 6 & CN1 != 6) {
          dt_tmp <- dt_tmp[(CN == CN1 & group_rho == gr1) | (CN == CN2 & group_rho == gr2), ]
          if (flag == "inter") {
            idxs <- dt_tmp[, .(
              s1 = sum(CN == CN1 & !duplicated(patient)),
              s2 = sum(CN == CN2 & !duplicated(patient))
            ), by = .(ID)]
          } else {
            idxs <- dt_tmp[, .(
              s1 = sum(CN == CN1 & (duplicated(patient) | duplicated(patient, fromLast = TRUE))),
              s2 = sum(CN == CN2 & (duplicated(patient) | duplicated(patient, fromLast = TRUE)))
            ), by = .(ID)]
          }
        }
        if (CN2 != 6 & CN1 == 6) {
          dt_tmp <- dt_tmp[(CN >= CN1 & group_rho == gr1) | (CN == CN2 & group_rho == gr2), ]
          if (flag == "inter") {
            idxs <- dt_tmp[, .(
              s1 = sum(CN >= CN1 & !duplicated(patient)),
              s2 = sum(CN == CN2 & !duplicated(patient))
            ), by = .(ID)]
          } else {
            idxs <- dt_tmp[, .(
              s1 = sum(CN >= CN1 & (duplicated(patient) | duplicated(patient, fromLast = TRUE))),
              s2 = sum(CN == CN2 & (duplicated(patient) | duplicated(patient, fromLast = TRUE)))
            ), by = .(ID)]
          }
        }
        if (CN2 == 6 & CN1 != 6) {
          dt_tmp <- dt_tmp[(CN == CN1 & group_rho == gr1) | (CN >= CN2 & group_rho == gr2), ]
          if (flag == "inter") {
            idxs <- dt_tmp[, .(
              s1 = sum(CN == CN1 & !duplicated(patient)),
              s2 = sum(CN >= CN2 & !duplicated(patient))
            ), by = .(ID)]
          } else {
            idxs <- dt_tmp[, .(
              s1 = sum(CN == CN1 & (duplicated(patient) | duplicated(patient, fromLast = TRUE))),
              s2 = sum(CN >= CN2 & (duplicated(patient) | duplicated(patient, fromLast = TRUE)))
            ), by = .(ID)]
          }
        }
        if (CN2 == 6 & CN1 == 6) {
          dt_tmp <- dt_tmp[(CN >= CN1 & group_rho == gr1) | (CN >= CN2 & group_rho == gr2), ]
          if (flag == "inter") {
            idxs <- dt_tmp[, .(
              s1 = sum(CN >= CN1 & !duplicated(patient)),
              s2 = sum(CN >= CN2 & !duplicated(patient))
            ), by = .(ID)]
          } else {
            idxs <- dt_tmp[, .(
              s1 = sum(CN >= CN1 & (duplicated(patient) | duplicated(patient, fromLast = TRUE))),
              s2 = sum(CN >= CN2 & (duplicated(patient) | duplicated(patient, fromLast = TRUE)))
            ), by = .(ID)]
          }
        }
      } else {
        dt_tmp <- dt_tmp[(group_rho == gr1) | (group_rho == gr2), ]
        if (flag == "inter") {
          idxs <- dt_tmp[, .(s1 = sum(!duplicated(patient))), by = .(ID)]
        } else {
          idxs <- dt_tmp[, .(s1 = sum(duplicated(patient) | duplicated(patient, fromLast = TRUE))), by = .(ID)]
        }
      }

      if (is.null(CN1)) {
        idxs <- idxs[s1 >= 2, (ID)]
      } else {
        if (CN2 == CN1) {
          idxs <- idxs[s1 >= 2, (ID)]
        } else {
          idxs <- idxs[s1 >= 1 & s2 >= 1, (ID)]
        }
      }

      dt_tmp <- dt_tmp[ID %in% idxs, ]
      dt_tmp[, ID2 := paste(ID, patient, sep = "_")]
      if (flag == "intra") {
        dt_tmp[, ID3 := paste(ID, LTXID, sep = "_")]
      }
      rm(idxs)

      if (!is.null(CN1)) {
        if (CN2 != 6 & CN1 != 6) {
          dt1 <- dt_tmp[CN == CN1 & group_rho == gr1, ]
          dt2 <- dt_tmp[CN == CN2 & group_rho == gr2, ]
        }
        if (CN2 != 6 & CN1 == 6) {
          dt1 <- dt_tmp[CN >= CN1 & group_rho == gr1, ]
          dt2 <- dt_tmp[CN == CN2 & group_rho == gr2, ]
        }
        if (CN2 == 6 & CN1 != 6) {
          dt1 <- dt_tmp[CN == CN1 & group_rho == gr1, ]
          dt2 <- dt_tmp[CN >= CN2 & group_rho == gr2, ]
        }
        if (CN2 == 6 & CN1 == 6) {
          dt1 <- dt_tmp[CN >= CN1 & group_rho == gr1, ]
          dt2 <- dt_tmp[CN >= CN2 & group_rho == gr2, ]
        }
      } else {
        dt1 <- dt_tmp[group_rho == gr1, ]
        dt2 <- dt_tmp[group_rho == gr2, ]
      }

      if (flag == "inter") {
        dt1 <- dt1[!duplicated(ID), ]
        dt1 <- dt1[!duplicated(ID2), ]
        dt2 <- dt2[!ID2 %in% dt1$ID2, ]
        dt2 <- dt2[!duplicated(ID)]
        dt2 <- dt2[ID %in% dt1$ID]
        dt1 <- dt1[ID %in% dt2$ID]
      } else {
        dt1 <- dt1[!duplicated(ID2), ]
        dt2 <- dt2[(!ID3 %in% dt1$ID3), ]
        dt2 <- dt2[!duplicated(ID2)]
        dt2 <- dt2[ID2 %in% dt1$ID2, ]
        dt1 <- dt1[ID2 %in% dt2$ID2, ]
        dt1[, ID := ID2]
        dt2[, ID := ID2]
      }

      cols <- c("ID", "m_t_low", "m_t_high", "m_t", "m_b_low", "m_b_high", "m_b")
      dt1[, c("CN1", "rho1", "id1", "gr1") := list(CN, rho, LTXID, group_rho)]
      dt2[, c("CN2", "rho2", "id2", "gr2") := list(CN, rho, LTXID, group_rho)]
      dt1 <- dt1[, .SD, .SDcols = c(cols, "CN1", "rho1", "id1", "gr1")]
      dt2 <- dt2[, .SD, .SDcols = c(cols, "CN2", "rho2", "id2", "gr2")]

      cols1 <- paste(colnames(dt1)[grepl("m_", colnames(dt1))], 1, sep = "_")
      cols2 <- paste(colnames(dt2)[grepl("m_", colnames(dt2))], 2, sep = "_")
      colnames(dt1)[grepl("m_", colnames(dt1))] <- cols1
      colnames(dt2)[grepl("m_", colnames(dt2))] <- cols2
      rm(cols1, cols2)

      df <- Reduce(
        function(...) merge(..., by = "ID", all = TRUE),
        list(dt1, dt2)
      )
      dt_tmp <- data.table(df)

      cat(nrow(dt_tmp), " CpGs to sample from\n")

      sample_size <- ifelse(!is.null(CN1), 1e4, 1e6)
      if (nrow(dt_tmp) >= 1e3) {
        if (nrow(dt_tmp) < 1e4) {
          dt_tmp <- dt_tmp[sample(seq_len(nrow(dt_tmp)), sample_size, replace = TRUE)]
        } else {
          dt_tmp <- dt_tmp[sample(seq_len(nrow(dt_tmp)), sample_size, replace = FALSE)]
        }
        rm(df)

        cols1 <- colnames(dt1)[grepl("m_t_[hl]", colnames(dt1))]
        cols2 <- colnames(dt2)[grepl("m_t_[hl]", colnames(dt2))]
        raw <- Overlap(
          dt_tmp[, .SD, .SDcols = cols1],
          dt_tmp[, .SD, .SDcols = cols2]
        ) == 0
        diffs <- dt_tmp[, m_t_1 - m_t_2]
        corr <- abs(diffs) >= effect_size & raw == TRUE
        rm(diffs)

        dt_tmp[, DMP_tt_raw := raw]
        dt_tmp[, DMP_tt_corr := corr]

        cols1 <- colnames(dt1)[grepl("m_b_[hl]", colnames(dt1))]
        cols2 <- colnames(dt2)[grepl("m_b_[hl]", colnames(dt2))]
        raw <- Overlap(
          dt_tmp[, .SD, .SDcols = cols1],
          dt_tmp[, .SD, .SDcols = cols2]
        ) == 0
        diffs <- dt_tmp[, m_b_1 - m_b_2]
        corr <- abs(diffs) >= effect_size & raw == TRUE
        rm(diffs)

        dt_tmp[, DMP_bb_raw := raw]
        dt_tmp[, DMP_bb_corr := corr]
        dt_list[[length(dt_list) + 1]] <- dt_tmp
        rm(raw, corr, cols1, cols2, dt1, dt2)
      }
    }
  }

  assign(paste0("dt_list_", flag), dt_list)
  cols_corrected <- grep("_corr", colnames(dt_results))
  rm(dt_list)

  dt_results <- rbindlist(dt_list)
  rm(dt_list)

  dt_results[, labels := ifelse(
    gr1 == 1 & gr2 == 1, "low_vs_low",
    ifelse(gr1 == 4 & gr2 == 4, "high_vs_high", "low_vs_high")
  )]

  labs <- unique(dt_results$labels)
  for (label in labs) {
    tmp <- dt_results[labels == label, ]
    dmp_cols <- colnames(tmp)[grepl("DMP_[tb][tb]_corr", colnames(tmp))]
    pdf_file <- file.path(
      out_dir,
      paste0(flag, "_", label, "_mt_mt_vs_mb_mb_DMPs_corrected_new.pdf")
    )
    pdf(pdf_file, width = 6, height = 6)
    print(plot(euler(tmp[, .SD, .SDcols = dmp_cols], shape = "ellipse"), quantities = TRUE))
    dev.off()
    Sys.sleep(1)
  }

  rm(dt_results, labs)
}

################################################################################
# Save list objects for downstream power analyses
################################################################################

save(dt_list_inter, dt_list_intra, file = rdata_out)
cat("Saved:", rdata_out, "\n")

# END
