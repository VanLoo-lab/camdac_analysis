#!/usr/bin/env Rscript

#' Figure 4I (DMP assignment bars) and 4J (clinical odds ratios) from cached cohort summaries.
#'
#' Expects \code{ASM_validation_Jul2021_v6.RData} (objects \code{df_summary3}, etc.)
#' produced by the Jul 2021 ASM plotting pipeline
#' (\code{R_scripts/Notes_ASM_plots_14Jul2021.R}, end of patient loop).

option_list <- list(
  optparse::make_option(
    "--validation_rdata", type = "character",
    default = "~/rdeocd/txMeth/ASM_validation_Jul2021_v6.RData"
  ),
  optparse::make_option(
    "--tx_dir", type = "character", default = "~/rdeocd/txMeth"
  ),
  optparse::make_option(
    "--work_dir", type = "character", default = "~/data/CAMDAC_paper/figure4"
  ),
  optparse::make_option(
    "--clinical_update", type = "character", default = "",
    help = paste0(
      "TRACERx methylation clinical table (pTNM stage); default ",
      "{tx_dir}/Tracerx_Methylation_Clinical_data_update.txt"
    )
  ),
  optparse::make_option("--panels", type = "character", default = "I,J"),
  optparse::make_option(
    "--clo_bi_contrast", type = "character", default = "other",
    help = paste0("Fig 4J third facet: 'other' = clonal bi-allelic vs all other DMPs ",
                  "(Notes_ASM_plots_14Jul2021.R ~L1149) or 'clo_mono' = vs clonal mono-allelic only")
  ),
  optparse::make_option(
    "--j_table", type = "character", default = "summary3",
    help = "Fig 4J cohort table: summary3 (default, Notes OR table) or summary2"
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

odds_ratio_2x2 <- function(mat) {
  if (requireNamespace("questionr", quietly = TRUE)) {
    as.numeric(questionr::odds.ratio(mat))
  } else {
    ft <- fisher.test(mat)
    c(as.numeric(ft$estimate), ft$conf.int[1], ft$conf.int[2], ft$p.value)
  }
}

pval_stars <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "")))
  )
}

validation_rdata <- path.expand(opt$validation_rdata)
if (!file.exists(validation_rdata)) {
  stop(
    "Missing cohort summary RData: ", validation_rdata, "\n",
    "Run Notes_ASM_plots_14Jul2021.R (summary_stats=TRUE) or provide --validation_rdata."
  )
}
load(validation_rdata)
if (!exists("df_summary3")) stop("df_summary3 not found in ", validation_rdata)

output_dir <- file.path(path.expand(opt$work_dir), "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panels <- toupper(trimws(strsplit(opt$panels, ",", fixed = TRUE)[[1]]))

tx_dir <- path.expand(opt$tx_dir)
clinical <- file.path(
  tx_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
clinical_update <- if (nzchar(opt$clinical_update)) {
  path.expand(opt$clinical_update)
} else {
  file.path(tx_dir, "Tracerx_Methylation_Clinical_data_update.txt")
}

merge_tracerx_clinical_stage <- function(df, clin_update_path) {
  if (!file.exists(clin_update_path)) {
    warning("Clinical update not found (no tumor stage for 4J): ", clin_update_path)
    df[, stage := NA_character_]
    return(df)
  }
  dt_clin <- as.data.table(read.delim(clin_update_path, stringsAsFactors = FALSE))
  dt_clin[, Patient := substr(Hospital_ID, 3, 9)]
  df[, pTNMStage_v8 := dt_clin[match(df$Patient, Patient), pTNMStage_v8]]
  df[, stage := ifelse(pTNMStage_v8 %in% c("1a", "1b"), "early", "late")]
  rec_col <- "Recurrence.1..yes.0.no"
  if (rec_col %in% names(dt_clin)) {
    df[, recurrence := dt_clin[match(df$Patient, Patient), get(rec_col)]]
  }
  df
}
if (file.exists(clinical)) {
  patient_info <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
  df_summary3[, hist := patient_info[match(df_summary3$Patient, Patient), Histology]]
  df_summary3[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
}

# Fig 4I colours: Notes_ASM_plots_14Jul2021.R ~L1327-1336 (fill alpha 0.5, same-colour outline).
# Named vectors (colours are matched by class name, not by factor order).
dmp_cols <- c(
  clo_AS = "deeppink3", clo_bi = "chartreuse4", sub_ambig = "deepskyblue2",
  sub_AS = "lightpink2", sub_bi = "gold2"
)
dmp_labels <- c(
  clo_AS = "clonal AS-DMP", clo_bi = "clonal bi-DMP", sub_ambig = "unclassified",
  sub_AS = "subclonal AS-DMP", sub_bi = "subclonal bi-DMP"
)

if ("I" %in% panels) {
  req_i <- c(
    "clo_AS_DMPs", "clo_bi_DMPs", "sub_AS_DMPs", "sub_bi_DMPs", "sub_ambig_DMPs"
  )
  if (!all(req_i %in% names(df_summary3))) {
    stop("df_summary3 missing DMP count columns for panel I: ",
         paste(setdiff(req_i, names(df_summary3)), collapse = ", "))
  }
  tmp <- data.table(
    Patient = df_summary3$Patient,
    pp = df_summary3$CRUKID,
    y = c(
      df_summary3$clo_AS_DMPs, df_summary3$clo_bi_DMPs, df_summary3$sub_AS_DMPs,
      df_summary3$sub_bi_DMPs, df_summary3$sub_ambig_DMPs
    ),
    coll = rep(c("clo_AS", "clo_bi", "sub_AS", "sub_bi", "sub_ambig"), each = nrow(df_summary3))
  )
  if (file.exists(clinical)) {
    pi <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
    tmp[, hist := pi[match(tmp$Patient, Patient), Histology]]
    tmp[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
  } else {
    tmp[, hist := "all"]
  }
  # order patients by clonal fraction of assigned DMPs (Notes ~L1319, ambiguous excluded)
  ordd <- df_summary3[, .(
    m = mean((clo_AS_DMPs + clo_bi_DMPs) /
      (clo_AS_DMPs + clo_bi_DMPs + sub_AS_DMPs + sub_bi_DMPs))
  ), by = .(CRUKID)]
  lv <- ordd$CRUKID[order(ordd$m)]
  tmp[, ID := factor(pp, levels = lv, labels = sub("^CRUK0*", "", lv), ordered = TRUE)]
  # stacking: first level on top (clonal AS top ... subclonal bi-allelic bottom)
  tmp[, coll := factor(coll, levels = names(dmp_cols))]

  theme_i <- theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 13, hjust = 0),
          axis.text.x = element_text(angle = 90, size = 8.5, vjust = 0.5, hjust = 1),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
          axis.title.y = element_text(size = 11), axis.text.y = element_text(size = 10),
          legend.position = "bottom", legend.text = element_text(size = 9),
          legend.key.size = unit(0.5, "lines"))
  scales_i <- list(
    scale_color_manual(name = "Epimutation class", values = dmp_cols, labels = dmp_labels, drop = FALSE),
    scale_fill_manual(name = "Epimutation class", values = dmp_cols, labels = dmp_labels, drop = FALSE)
  )
  # top: counts incl. unassigned; bottom: relative ratio excluding unassigned (sub_ambig)
  p_counts <- ggplot(tmp, aes(x = ID, y = y, col = coll, fill = coll)) +
    facet_grid(. ~ hist, scales = "free_x", space = "free_x") + xlab("") +
    geom_col(alpha = 0.5, linewidth = 0.2) + scales_i +
    scale_y_continuous(name = "DMP assignments\n(inc. unassigned)") + theme_i
  p_ratio <- ggplot(tmp[coll != "sub_ambig"], aes(x = ID, y = y, col = coll, fill = coll)) +
    facet_grid(. ~ hist, scales = "free_x", space = "free_x") + xlab("") +
    geom_col(alpha = 0.5, linewidth = 0.2, position = "fill") + scales_i +
    scale_y_continuous(name = "relative ratio\n(exc. unassigned)") +
    theme_i + theme(strip.text = element_blank())
  leg_guides <- guides(col = guide_legend(nrow = 2, byrow = TRUE, title.position = "left"),
                       fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "left"))
  # legend taken from the counts plot so that "unclassified" (absent from the ratio panel) is shown
  legend_i <- cowplot::get_legend(p_counts + leg_guides)
  p_ratio_leg <- p_ratio + leg_guides
  p_counts_clean <- p_counts + theme(legend.position = "none", axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank())
  p_ratio_clean <- p_ratio + theme(legend.position = "none")
  fig_i <- cowplot::plot_grid(p_counts_clean, p_ratio_clean, legend_i, ncol = 1, align = "v",
                              axis = "lr", rel_heights = c(1, 1.3, 0.28))
  ggsave(file.path(output_dir, "Fig4I_DMP_assignments.pdf"), fig_i, width = 8.2, height = 6.2)
  ggsave(file.path(output_dir, "Fig4I_DMP_assignments.png"), fig_i, width = 8.2, height = 6.2, dpi = 160)
  ggsave(file.path(output_dir, "Fig4I_DMP_assignments_counts.pdf"), p_counts, width = 8.2, height = 3.4)
  ggsave(file.path(output_dir, "Fig4I_DMP_assignments_ratio.pdf"), p_ratio_leg, width = 8.2, height = 3.6)
  fwrite(tmp, file.path(output_dir, "Fig4I_DMP_assignments.tsv"), sep = "\t")
}

if ("J" %in% panels) {
  j_table <- match.arg(opt$j_table, c("summary3", "summary2"))
  if (j_table == "summary2" && exists("df_summary2")) {
    df <- copy(df_summary2)
    message("Fig 4J: df_summary2 (per-region rows; n=", nrow(df), ")")
  } else {
    df <- copy(df_summary3)
    message(
      "Fig 4J: df_summary3 (one row per patient; Notes OR input; n=", nrow(df), ")"
    )
  }
  if (!file.exists(clinical)) {
    stop("Panel J requires clinical table: ", clinical)
  }
  if (!exists("patient_info")) {
    patient_info <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
  }
  patient_info <- patient_info[!grepl("LN", Sample)]
  df[, hist := patient_info[match(df$Patient, Patient), Histology]]
  df[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
  wgd_map <- patient_info[, .(WGD_n = round(mean(as.numeric(WGD_reseg)), 0)), by = Patient]
  df[, WGD := wgd_map[match(df$Patient, Patient), WGD_n]]
  df[, WGD := ifelse(WGD == 0, "non-WGD", "WGD")]
  df <- merge_tracerx_clinical_stage(df, clinical_update)
  ttmp <- patient_info[, .(num_sample = length(Sample)), by = Patient]
  df[, num_sample := ttmp[match(df$Patient, Patient), num_sample]]
  # Paper Fig. 4J: high > 3 regions sequenced, low otherwise (not cohort median).
  df[, num_samples := ifelse(num_sample > 3, "high", "low")]
  sbs_file <- file.path(
    tx_dir, "tx421_SBS_signatures", "tx421_signature_weights_perTumour.txt"
  )
  if (file.exists(sbs_file)) {
    dt_SBS <- as.data.table(read.delim(sbs_file, stringsAsFactors = FALSE))
    dt_SBS[, Patient := substr(TumourID, 1, 6)]
    df[, SBS4 := dt_SBS[match(df$Patient, Patient), SBS4]]
    df[, SBS4 := ifelse(SBS4 < 0.35, "low", "high")]
  }
  n <- nrow(df)
  # Notes OR (~L1177): mean(tmp_col2) by clinical level × DMP-type pair.
  # Long table: two rows per patient (or per region if summary2). Stratum
  # cells are therefore mean DMP count **per patient** (summary3 default),
  # not cohort sums and not conflated with mean across regions unless using
  # summary2 (see README — cached summary2 may lack subclonal counts).
  tmp <- data.table(
    Patient = rep(df$Patient, 2L),
    hist = rep(df$hist, 2L),
    stage = rep(df$stage, 2L),
    recurrence = rep(df$recurrence, 2L),
    SBS4 = rep(df$SBS4, 2L),
    WGD = rep(df$WGD, 2L),
    num_samples = rep(df$num_samples, 2L),
    clonality = c(df$clo_AS_DMPs + df$clo_bi_DMPs, df$sub_AS_DMPs + df$sub_bi_DMPs),
    coll1 = rep(c("clonal", "subclonal"), each = n),
    allele = c(df$clo_AS_DMPs + df$sub_AS_DMPs, df$clo_bi_DMPs + df$sub_bi_DMPs),
    coll2 = rep(c("mono", "bi"), each = n),
    clo_bi = c(df$clo_bi_DMPs, if (opt$clo_bi_contrast == "other") {
      df$sub_AS_DMPs + df$sub_bi_DMPs + df$clo_AS_DMPs
    } else df$clo_AS_DMPs),
    coll4 = rep(c("clo_bi", "other"), each = n)
  )
  tmp[, coll1 := factor(coll1, levels = c("clonal", "subclonal"), ordered = TRUE)]
  tmp[, coll4 := factor(coll4, levels = c("clo_bi", "other"), ordered = TRUE)]

  # Covariates: histology, stage, recurrence, WGD, SBS4, number of samples.
  covars <- c("hist", "stage", "recurrence", "SBS4", "WGD", "num_samples")
  dt_odds_list <- list()
  for (yy in covars) {
    for (flagg in c("allele", "clonality", "clo_bi")) {
      lab <- ifelse(flagg == "allele", "coll2",
                    ifelse(flagg == "clonality", "coll1", "coll4"))
      dt_odds <- data.table(y = yy, flag = flagg)
      tmp[, tmp_col1 := get(yy)]
      tmp[, tmp_col2 := get(flagg)]
      tmp[, tmp_col3 := get(lab)]
      mat_tmp <- tmp[!is.na(tmp_col1), .(
        x = mean(tmp_col2),
        n_units = .N,
        n_patients = uniqueN(Patient)
      ), by = .(tmp_col3, tmp_col1)]
      if (nrow(mat_tmp) < 4L) next
      mat <- matrix(mat_tmp$x, nrow = 2)
      if (nrow(mat) != 2 || ncol(mat) != 2) next
      # Paper Fig. 4J: histology OR is LUAD vs LUSC (swap row order vs default sort).
      if (yy == "hist") {
        mat <- mat[rev(seq_len(nrow(mat))), , drop = FALSE]
      }
      tmp_od <- odds_ratio_2x2(mat)
      dt_odds[, OR := tmp_od[1]]
      dt_odds[, OR_min := tmp_od[2]]
      dt_odds[, OR_max := tmp_od[3]]
      dt_odds[, pval := tmp_od[4]]
      dt_odds[, j_table := j_table]
      dt_odds[, n_units := mat_tmp[, sum(n_units)]]
      dt_odds[, n_patients := uniqueN(tmp[!is.na(tmp_col1), Patient])]
      dt_odds_list[[length(dt_odds_list) + 1L]] <- dt_odds
    }
  }
  dt_odds <- rbindlist(dt_odds_list, fill = TRUE)
  dt_odds[, contrast := data.table::fcase(
    flag == "allele", "mono-allelic vs bi-allelic DMPs",
    flag == "clonality", "clonal vs subclonal DMPs",
    flag == "clo_bi", "clonal bi-allelic vs clonal mono-allelic DMPs"
  )]
  dt_odds[, stars := pval_stars(pval)]
  dt_odds[, flag_label := factor(
    flag,
    levels = c("allele", "clonality", "clo_bi"),
    labels = c("Allele-specific", "Clonal", "Clonal bi-allelic")
  )]
  dt_odds[, y_label := factor(
    y,
    levels = covars,
    labels = c(
      "histology",
      "tumor stage",
      "recurrence",
      "SBS4",
      "WGD",
      "number of\nsamples"
    )
  )]
  or_x_max <- max(
    2,
    ceiling(max(dt_odds$OR_max, dt_odds$OR, na.rm = TRUE) * 1.08 * 10) / 10
  )
  dt_odds[, star_x := pmin(pmax(OR_max, OR) * 1.04, or_x_max * 0.98)]
  pval_key <- data.table(
    flag_label = factor("Clonal bi-allelic", levels = levels(dt_odds$flag_label)),
    y_label = factor("number of\nsamples", levels = levels(dt_odds$y_label)),
    x = 1.05,
    label = "* 0.05\n** 0.01\n*** 0.001"
  )
  p <- ggplot(dt_odds) +
    theme_classic(base_size = 9) +
    geom_col(
      aes(y = y_label, x = OR),
      width = 0.72,
      fill = "grey82",
      colour = "grey20",
      linewidth = 0.25
    ) +
    geom_linerange(
      aes(y = y_label, xmin = OR_min, xmax = OR_max),
      linewidth = 0.35,
      colour = "grey20"
    ) +
    geom_text(
      aes(y = y_label, x = star_x, label = stars),
      data = dt_odds[stars != ""],
      size = 3,
      hjust = 0,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = pval_key,
      aes(x = x, y = y_label, label = label),
      inherit.aes = FALSE,
      size = 2.5,
      hjust = 0,
      vjust = 0,
      colour = "grey15",
      lineheight = 0.9
    ) +
    facet_wrap(~flag_label, nrow = 1, scales = "fixed") +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = "grey40") +
    scale_x_continuous(
      "odds ratio",
      limits = c(0, or_x_max),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_discrete(name = NULL, limits = rev(levels(dt_odds$y_label))) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 9),
      panel.background = element_blank(),
      axis.text.y = element_text(size = 8, hjust = 1)
    )
  ggsave(file.path(output_dir, "Fig4J_clinical_odds_ratios.pdf"), p,
         width = 7.5, height = 3.2)
  fwrite(dt_odds, file.path(output_dir, "Fig4J_clinical_odds_ratios.tsv"), sep = "\t")

  # Published Fig 4J values read off the paper figure (approximate, +/- ~0.05).
  # order within each facet = covars: hist, stage, recurrence, SBS4, WGD, num_samples
  ref <- data.table(
    flag = rep(c("allele", "clonality", "clo_bi"), each = 6L),
    y = rep(c("hist", "stage", "recurrence", "SBS4", "WGD", "num_samples"), 3L),
    OR_published = c(
      0.68, 1.60, 0.85, 1.00, 1.03, 1.72,
      1.15, 1.00, 0.93, 1.20, 1.33, 1.30,
      1.80, 0.75, 0.72, 0.75, 0.55, 0.92
    ),
    stars_published = c(
      "***", "***", "*", "", "*", "***",
      "***", "", "***", "***", "***", "***",
      "***", "***", "***", "***", "***", ""
    )
  )
  qc <- merge(ref, dt_odds[, .(flag, y, OR, OR_min, OR_max, pval, stars)], by = c("flag", "y"), all = TRUE)
  qc[, flag := factor(flag, levels = c("allele", "clonality", "clo_bi"))]
  qc[, y := factor(y, levels = covars)]
  setorder(qc, flag, y)
  qc[, OR_delta := OR - OR_published]
  qc[, within_0.15 := abs(OR_delta) <= 0.15]
  qc[, same_direction := (OR > 1) == (OR_published > 1) | abs(OR_published - 1) < 0.06]
  qc[, same_stars := stars == stars_published]
  fwrite(qc, file.path(output_dir, "Fig4J_qc_vs_published_panel.tsv"), sep = "\t")
}

cat("Figure 4 cohort panels written to", output_dir, "\n")
