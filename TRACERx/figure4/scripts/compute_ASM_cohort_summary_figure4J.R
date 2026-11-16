#!/usr/bin/env Rscript
#
# Refresh ASM cohort summary tables for Figure 4I/J without overwriting
# ASM_validation_Jul2021_v6.RData.
#
# Matches Notes_ASM_plots_14Jul2021.R sample loop (flag3, hyper methylation,
# multip="all" per sample) but uses CpG-level flag3 counts (not UMAP cluster
# class). df_summary3 = one row per patient (sum of sample-level counts).

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
  optparse::make_option(
    "--tx_dir", type = "character", default = "~/rdeocd/txMeth"
  ),
  optparse::make_option(
    "--out_rdata", type = "character",
    default = "~/data/CAMDAC_paper/figure4/outputs/ASM_cohort_summary_figure4J_refresh.RData"
  ),
  optparse::make_option("--meth", type = "character", default = "hyper")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

suppressPackageStartupMessages({
  library(data.table)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(file_arg))
} else {
  path.expand("~/data/CAMDAC_paper/figure4/scripts")
}
source(file.path(script_dir, "asm_utils.R"))

tx_dir <- path.expand(opt$tx_dir)
meth <- opt$meth
thres <- ifelse(meth == "hyper", 0.1, 0.8)

clinical <- file.path(
  tx_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
if (!file.exists(clinical)) {
  stop("Clinical table not found: ", clinical)
}
patient_info <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
patient_info <- patient_info[!grepl("LN", Sample)]
patient_ids <- unique(patient_info$Patient)
patient_ids <- patient_ids[patient_ids != resolve_tx_patient_id("CRUK0047", .tx_map)]

assign_flag3 <- function(dt) {
  dt[, is_DM_t12 := ((m_t_1_low > m_t_2_high) | (m_t_2_low > m_t_1_high)) &
       abs(m_diff_t12) > 0.25]
  dt[, flag3 := ifelse(
    !is.na(pwr_tt) & pwr_tt == TRUE & is_DM_t12 == TRUE & m_t_1 > m_t_2,
    "clonal-AS-DMP1",
    ifelse(
      !is.na(pwr_tt) & pwr_tt == TRUE & is_DM_t12 == TRUE & m_t_1 < m_t_2,
      "clonal-AS-DMP2",
      ifelse(
        !is.na(pwr_tn) & pwr_tn == TRUE & !is.na(DMP_t) & abs(m_diff_t12) <= 0.25 &
          abs((m_t_1 + m_t_2) / 2 - m_n) > 0.65,
        "clonal-bi-DMP",
        ifelse(
          !is.na(pwr_tn) & pwr_tn == TRUE & !is.na(DMP_t) &
            !((m_t_1_low > m_t_2_high) | (m_t_2_low > m_t_1_high)),
          "subclonal-bi-DMP",
          ifelse(
            !is.na(pwr_tt) & pwr_tt == TRUE &
              ((m_t_1_low > m_t_2_high) | (m_t_2_low > m_t_1_high)) &
              m_diff_t12 > abs((m_t_1 + m_t_2) / 2 - m_n) & m_t_1 > m_t_2,
            "subclonal-AS-DMP1",
            ifelse(
              !is.na(pwr_tt) & pwr_tt == TRUE &
                ((m_t_1_low > m_t_2_high) | (m_t_2_low > m_t_1_high)) &
                m_diff_t12 > abs((m_t_1 + m_t_2) / 2 - m_n) & m_t_1 < m_t_2,
              "subclonal-AS-DMP2",
              ifelse(
                !is.na(pwr_tt) & !is.na(pwr_tn) & (pwr_tt == FALSE | pwr_tn == FALSE),
                "no-power",
                "non-DMP"
              )
            )
          )
        )
      )
    )
  )]
  dt[, flag3 := factor(
    flag3,
    levels = c(
      "clonal-AS-DMP1", "clonal-AS-DMP2", "clonal-bi-DMP",
      "subclonal-AS-DMP1", "subclonal-AS-DMP2", "subclonal-bi-DMP",
      "non-DMP", "no-power"
    ),
    ordered = TRUE
  )]
  dt
}

sample_summary_row <- function(patient_id, sample_id, dt_plot) {
  data.table(
    Patient = patient_id,
    CRUKID = patient_info[Patient == patient_id, CRUKID[1]],
    Sample = sample_id,
    multip = "all",
    flag = meth,
    WGD = patient_info[Patient == patient_id & Sample == sample_id, WGD_reseg],
    no_power = dt_plot[, sum(flag3 == "no-power")],
    no_power_clo_AS = dt_plot[, sum(flag3 == "no-power" & abs(m_diff_t12) > 0.4)],
    no_power_sub_bi_AS = dt_plot[, sum(
      flag3 == "no-power" &
        abs(m_diff_t12) <= 0.4 & abs(m_t_corr - m_n) > 0.2
    )],
    clo_bi_DMPs = dt_plot[, sum(flag3 == "clonal-bi-DMP")],
    clo_AS_DMPs = dt_plot[, sum(flag3 %in% c("clonal-AS-DMP1", "clonal-AS-DMP2"))],
    sub_bi_DMPs = dt_plot[, sum(flag3 == "subclonal-bi-DMP")],
    sub_AS_DMPs = dt_plot[, sum(flag3 %in% c("subclonal-AS-DMP1", "subclonal-AS-DMP2"))],
    sub_ambig_DMPs = 0L
  )
}

add_perc_cols <- function(dt) {
  dt[, perc_IM_AS_v_sub := (clo_AS_DMPs + sub_AS_DMPs) /
       (clo_AS_DMPs + sub_bi_DMPs + sub_AS_DMPs)]
  dt[, perc_AS_v_bi := (clo_AS_DMPs + sub_AS_DMPs) /
       (clo_AS_DMPs + sub_bi_DMPs + sub_AS_DMPs + clo_bi_DMPs)]
  dt[, perc_clo_v_sub := (clo_AS_DMPs + clo_bi_DMPs) /
       (clo_AS_DMPs + clo_bi_DMPs + sub_bi_DMPs + sub_AS_DMPs)]
  invisible(dt)
}

df_list3 <- list()
skipped <- list()

for (patient_id in patient_ids) {
  sample_ids <- patient_info[Patient == patient_id, unique(Sample)]
  sample_ids <- sample_ids[order(
    patient_info[Patient == patient_id, cov_t],
    decreasing = TRUE
  )]
  for (sample_id in sample_ids) {
    rdata <- file.path(
      tx_dir, patient_id, "Allele_specific_analysis", sample_id,
      paste0("DMP_classification_", sample_id, ".RData")
    )
    if (!file.exists(rdata)) {
      skipped[[length(skipped) + 1L]] <- data.table(
        Patient = patient_id, Sample = sample_id, reason = "missing_RData"
      )
      next
    }
    load(rdata)
    dt <- dt_DMP_classification
    rm(dt_DMP_classification)
    dt <- assign_flag3(dt)
    cols2 <- c(
      "chrom", "start", "end", "m_n", "m_n_low", "m_n_high",
      "m_b", "m_1", "m_2", "m_t_corr", "m_t_1_raw", "m_t_2_raw",
      "m_t_1", "m_t_2", "m_diff_t12", "n1", "n2", "multip", "flag", "flag3"
    )
    ov <- dt[n1 != 0 & n2 != 0, .SD, .SDcols = cols2]
    rm(dt)
    if (thres < 0.5) {
      dt_plot <- ov[m_n_high < thres]
    } else {
      dt_plot <- ov[m_n_low > thres]
    }
    rm(ov)
    if (nrow(dt_plot) == 0L) {
      skipped[[length(skipped) + 1L]] <- data.table(
        Patient = patient_id, Sample = sample_id, reason = "empty_dt_plot"
      )
      next
    }
    row <- sample_summary_row(patient_id, sample_id, dt_plot)
    add_perc_cols(row)
    df_list3[[length(df_list3) + 1L]] <- row
    message(patient_id, "-", sample_id, " CpGs=", nrow(dt_plot))
    rm(dt_plot, row)
  }
}

if (length(df_list3) == 0L) stop("No sample summaries produced.")
df_summary2 <- rbindlist(df_list3, use.names = TRUE, fill = TRUE)
rm(df_list3)

count_cols <- c(
  "no_power", "no_power_clo_AS", "no_power_sub_bi_AS",
  "clo_bi_DMPs", "clo_AS_DMPs", "sub_bi_DMPs", "sub_AS_DMPs", "sub_ambig_DMPs"
)
df_summary3 <- df_summary2[, lapply(.SD, sum), .SDcols = count_cols, by = .(Patient, CRUKID, flag)]
wgd_patient <- patient_info[, .(
  WGD = round2(mean(as.numeric(WGD_reseg)), digits = 0)
), by = Patient]
df_summary3[, WGD := wgd_patient[match(df_summary3$Patient, wgd_patient$Patient), WGD]]
add_perc_cols(df_summary3)

refresh_meta <- list(
  built = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  source_script = "compute_ASM_cohort_summary_figure4J.R",
  notes_reference = "Notes_ASM_plots_14Jul2021.R (flag3, hyper, multip=all)",
  count_definition = paste(
    "CpG-level counts from flag3 on dt_plot;",
    "df_summary2 = one row per tumour sample;",
    "df_summary3 = sum of count columns across samples per patient;",
    "sub_ambig_DMPs = 0 (UMAP-only in original Notes df_summary3)"
  ),
  n_patients_clinical = length(patient_ids),
  n_patients_summary3 = nrow(df_summary3),
  n_rows_summary2 = nrow(df_summary2),
  meth = meth
)
skipped_samples <- if (length(skipped)) rbindlist(skipped) else data.table()

out_rdata <- path.expand(opt$out_rdata)
protected <- path.expand("~/rdeocd/txMeth/ASM_validation_Jul2021_v6.RData")
if (normalizePath(out_rdata, mustWork = FALSE) ==
    normalizePath(protected, mustWork = FALSE)) {
  stop("Refusing to write to protected v6 RData: ", protected)
}
dir.create(dirname(out_rdata), recursive = TRUE, showWarnings = FALSE)
save(
  df_summary2, df_summary3, refresh_meta, skipped_samples,
  file = out_rdata
)
message("Wrote ", out_rdata)
message(
  "df_summary2: ", nrow(df_summary2), " rows; df_summary3: ",
  nrow(df_summary3), " patients"
)
