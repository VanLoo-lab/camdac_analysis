#!/usr/bin/env Rscript

#' Figure 4A — per-patient pie charts of imprinting state at known loci.
#'
#' Source logic: \code{R_scripts/Notes_imprinted_genes_v2.R} (~L729–1025),
#' data: \code{imprinting_summary_tumours_Oct2021_v3.RData} (\code{dt_IM_tumours}).

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
  optparse::make_option(
    "--tumours_rdata", type = "character", default = "",
    help = "Override path to imprinting_summary_tumours_Oct2021_v3.RData"
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

tx_dir <- path.expand(opt$tx_dir)
work_dir <- path.expand(opt$work_dir)
output_dir <- file.path(work_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumours_rdata <- if (nzchar(opt$tumours_rdata)) {
  path.expand(opt$tumours_rdata)
} else {
  file.path(tx_dir, "imprinting_summary_tumours_Oct2021_v3.RData")
}
if (!file.exists(tumours_rdata)) stop("Missing: ", tumours_rdata)

clinical <- file.path(
  tx_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
if (!file.exists(clinical)) stop("Missing clinical table: ", clinical)
patient_info <- as.data.table(read.delim(clinical, stringsAsFactors = FALSE))
patient_info <- patient_info[!grepl("LN", Sample)]
patient_info[, LTX_ID := paste0(Patient, "_", Sample)]

load(tumours_rdata)
dt_IM_tumours <- as.data.table(dt_IM_tumours)

select_imprinting_bin_ids <- function(dt) {
  tmp <- dt[
    Status == "Imprinted" & flag != "no_info",
    .(most_abberated_bin_patient = sum(!grepl("Imprinted", flag))),
    by = .(id, Gene, Status, patient, chrom)
  ]
  tmp[, max_number_abberated := max(most_abberated_bin_patient),
      by = .(Gene, Status, patient, chrom)]
  tmp[, sum_abberated := sum(most_abberated_bin_patient),
      by = .(id, Gene, Status, chrom)]
  tmp <- tmp[most_abberated_bin_patient == max_number_abberated]
  tmp2 <- tmp[, .(
    ids = ifelse(length(id) == 1L, id, id[which.max(sum_abberated)])
  ), by = .(Gene, Status, patient, chrom)]
  tmp2[, paste(patient, ids, sep = "_")]
}

build_dt_IM_tumours_ps <- function(dt, ids) {
  sub <- dt[paste(patient, id, sep = "_") %in% ids]
  ps <- sub[Status == "Imprinted", .(
    id = paste(sort(unique(id)), collapse = ","),
    driver = ifelse(sum(as.numeric(driver)) > 0, "y", "n"),
    n1 = mean(n1, na.rm = TRUE),
    n2 = mean(n2, na.rm = TRUE),
    m_diff_n12 = mean(m_diff_n12, na.rm = TRUE),
    m_diff_t12 = mean(m_diff_t12, na.rm = TRUE),
    m_t_1 = mean(m_t_1, na.rm = TRUE),
    m_t_2 = mean(m_t_2, na.rm = TRUE),
    m_n_1 = mean(m_n_1, na.rm = TRUE),
    m_n_2 = mean(m_n_2, na.rm = TRUE),
    CpGs_n = ifelse(sum(!is.na(CpGs_n)) > 0, mean(CpGs_n, na.rm = TRUE), NA_real_),
    CpGs_t = ifelse(sum(!is.na(CpGs_t)) > 0, mean(CpGs_t, na.rm = TRUE), NA_real_),
    m_t = ifelse(
      sum(flag == "LOH") > 0, mean(m_t[grepl("LOH", flag)], na.rm = TRUE),
      ifelse(
        sum(grepl("Imprinted", flag)) == .N |
          (sum(flag == "no_info") > 0 &
             sum(grepl("Imprinted", flag)) == sum(flag != "no_info")),
        NA_real_,
        ifelse(sum(grepl("LOI", flag)) > 0,
               mean(m_t[grepl("LOI", flag)], na.rm = TRUE),
               mean(m_t, na.rm = TRUE))
      )
    ),
    mean_m_n_phased = mean((m_n_1 + m_n_2) / 2, na.rm = TRUE),
    mean_m_t_phased = ifelse(
      sum(flag == "LOH") > 0,
      mean(c(
        m_t_1[grepl("LOH", flag) & !is.na(m_t_1) & n1 != 0],
        m_t_2[grepl("LOH", flag) & !is.na(m_t_2) & n2 != 0]
      ), na.rm = TRUE),
      ifelse(
        sum(grepl("Imprinted", flag)) == .N |
          (sum(flag == "no_info") > 0 &
             sum(grepl("Imprinted", flag)) == sum(flag != "no_info")),
        NA_real_,
        ifelse(
          sum(grepl("LOI", flag)) > 0,
          mean(c((m_t_1 + m_t_2) / 2)[grepl("LOI", flag)], na.rm = TRUE),
          mean(c((m_t_1 + m_t_2) / 2), na.rm = TRUE)
        )
      )
    ),
    flag = ifelse(
      sum(flag == "LOH") > 0, "LOH",
      ifelse(
        sum(grepl("Imprinted", flag)) == .N, "Imprinted",
        ifelse(
          sum(grepl("LOI", flag)) == .N, "LOI",
          ifelse(
            sum(grepl("LOI", flag)) > 0, "partial_LOI",
            ifelse(sum(flag == "no_info") > 0, "no_info", "investigate")
          )
        )
      )
    )
  ), by = .(Gene, Status, patient, sample, chrom)]

  ps[, M_m_t_phased2 := ifelse(
    grepl("Imprinted", flag), NA_real_,
    ifelse(mean_m_t_phased > mean_m_n_phased, 1, 0)
  )]
  ps[, LTX_ID := paste0(patient, "_", sample)]
  ps[, CRUKID := patient_info$CRUKID[match(ps$LTX_ID, patient_info$LTX_ID)]]
  ps[is.na(CRUKID), CRUKID := patient]
  ps
}

apply_manual_imprinting_fixes <- function(ps) {
  ps[patient == resolve_tx_patient_id("CRUK0071", .tx_map) & grepl("GNAS", Gene) & sample == "R1", `:=`(flag = "LOI", M_m_t_phased2 = 1)]
  ps[patient == resolve_tx_patient_id("CRUK0046", .tx_map) & grepl("MIMT1", Gene) & sample == "R4", `:=`(flag = "LOI", M_m_t_phased2 = 0)]
  ps[patient == resolve_tx_patient_id("CRUK0083", .tx_map) & grepl("H19", Gene) & sample == "R2", `:=`(flag = "LOI", M_m_t_phased2 = 0)]
  ps[patient == resolve_tx_patient_id("CRUK0049", .tx_map) & grepl("H19", Gene) & sample == "R2", `:=`(flag = "Imprinted", M_m_t_phased2 = NA_real_)]
  ps[patient == resolve_tx_patient_id("CRUK0049", .tx_map) & grepl("SNRPN", Gene) & sample == "R1", `:=`(flag = "Imprinted", M_m_t_phased2 = NA_real_)]
  ps
}

# Published Fig. 4A palette (solid grey Imprinted; LOI/LOH UM/M colors).
imprinting_levels <- c("Imprinted", "LOI_UM", "LOI_M", "LOH_UM", "LOH_M")
imprinting_fill <- c(
  Imprinted = "grey82",
  LOI_UM = "darkolivegreen3",
  LOI_M = "forestgreen",
  LOH_UM = "steelblue3",
  LOH_M = "hotpink3"
)
imprinting_labels <- c(
  Imprinted = "Imprinted",
  LOI_UM = "LOI, UM",
  LOI_M = "LOI, M",
  LOH_UM = "LOH, UM",
  LOH_M = "LOH, M"
)

pie_tables_for_gene <- function(ps, gene_pattern) {
  d <- ps[grepl(gene_pattern, Gene)]
  if (!nrow(d)) return(NULL)
  tmp2 <- d[, .(
    all = length(unique(LTX_ID)),
    Imprinted = length(unique(LTX_ID[grepl("Imprinted", flag)])) / length(unique(LTX_ID)),
    LOI_UM = length(unique(LTX_ID[grepl("LOI", flag) & M_m_t_phased2 == 0])) / length(unique(LTX_ID)),
    LOI_M = length(unique(LTX_ID[grepl("LOI", flag) & M_m_t_phased2 == 1])) / length(unique(LTX_ID)),
    LOH_UM = length(unique(LTX_ID[grepl("LOH", flag) & M_m_t_phased2 == 0])) / length(unique(LTX_ID)),
    LOH_M = length(unique(LTX_ID[grepl("LOH", flag) & M_m_t_phased2 == 1])) / length(unique(LTX_ID))
  ), by = .(Gene, CRUKID)]
  tmp2 <- tmp2[all > 1]
  if (!nrow(tmp2)) return(NULL)

  tmp3 <- melt(
    tmp2,
    id.vars = c("CRUKID", "Gene"),
    measure.vars = imprinting_levels
  )
  tmp3[, variable := factor(variable, levels = imprinting_levels)]
  tmp3[, all := tmp2$all[match(paste(CRUKID, Gene), paste(tmp2$CRUKID, tmp2$Gene))]]
  tmp3[, patient_num := as.integer(gsub("^CRUK0*", "", CRUKID))]
  ord_levels <- tmp3[
    variable == "Imprinted",
    unique(patient_num[order(value, all)])
  ]
  tmp3[, patient_ord := factor(patient_num, levels = ord_levels, ordered = TRUE)]
  tmp3 <- tmp3[order(patient_ord)]
  slices <- tmp2[, .(pie_slices = seq(0, 1, 1 / all)), by = .(CRUKID, Gene)]
  slices[, patient_ord := tmp3$patient_ord[match(
    paste(CRUKID, Gene),
    paste(tmp3$CRUKID, tmp3$Gene)
  )]]
  slices <- unique(slices[!is.na(patient_ord)])
  list(tmp3 = tmp3, slices = slices)
}

plot_gene_pies <- function(tables, title, ncol = 6) {
  tmp3 <- tables$tmp3
  slices <- tables$slices
  ggplot(tmp3, aes(x = "", y = value, fill = variable)) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 6),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 9, hjust = 0)
    ) +
    geom_bar(width = 1, stat = "identity", colour = NA) +
    geom_hline(
      data = slices,
      aes(yintercept = pie_slices, group = patient_ord),
      col = "black",
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = imprinting_fill,
      breaks = imprinting_levels,
      labels = imprinting_labels,
      name = NULL
    ) +
    coord_polar("y", start = 0) +
    facet_wrap(~patient_ord, ncol = ncol) +
    ggtitle(title)
}

ids <- select_imprinting_bin_ids(dt_IM_tumours)
dt_ps <- apply_manual_imprinting_fixes(build_dt_IM_tumours_ps(dt_IM_tumours, ids))

gene_specs <- list(
  list(pat = "H19", title = "IGF2/H19", ncol = 11),
  list(pat = "GNAS", title = "GNAS", ncol = 5),
  list(pat = "MIMT1,PEG3,ZIM2", title = "MIMT1/PEG3/ZIM2", ncol = 4),
  list(pat = "BLCAP,NNAT", title = "BLCAP/NNAT", ncol = 4),
  list(pat = "KCNQ1", title = "KCNQ1", ncol = 4),
  list(pat = "SNRPN", title = "SNRPN", ncol = 2)
)

plots <- list()
for (spec in gene_specs) {
  tab <- pie_tables_for_gene(dt_ps, spec$pat)
  if (is.null(tab)) next
  plots[[spec$title]] <- plot_gene_pies(tab, spec$title, ncol = spec$ncol)
}

if (!length(plots)) stop("No imprinting pie panels produced.")

out_pdf <- file.path(output_dir, "Fig4A_imprinting_overview.pdf")
ggsave(
  out_pdf,
  gridExtra::arrangeGrob(grobs = plots, ncol = 1),
  width = 14,
  height = 16,
  limitsize = FALSE
)
message("Wrote ", out_pdf)

for (nm in names(plots)) {
  ggsave(
    file.path(output_dir, paste0("Fig4A_", gsub("/", "_", nm), "_pies.pdf")),
    plots[[nm]],
    width = 8,
    height = 3.5,
    limitsize = FALSE
  )
}
