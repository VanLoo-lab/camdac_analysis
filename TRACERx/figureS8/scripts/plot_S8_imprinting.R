#' @title Supplementary Figure S8C (NSCLC imprinting overview)
#' @author Lili Larose Cadieux
#' @description Cohort overview of somatic changes at imprinted loci
#' (LOI / LOH / Imprinted) from TRACERx RRBS. Labels use CRUKID.
#'
#' Source: \code{Notes_imprinted_genes_v2.R}; cached
#' \code{imprinting_summary_gene_subset_for_plot_Oct2021.RData}.
#'
#' S8A–B (locus tracks for IGF2/H19 and MEST) and S8D (CRC WGBS overview)
#' need allele-specific methylation tracks / CRC inputs — see README.

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

resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) return(path)
  file.path(base_dir, path)
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

load_ltx_cruk_map <- function(clinical_info) {
  dt <- unique(fread(path.expand(clinical_info))[, .(Patient, CRUKID)])
  dt[!is.na(Patient) & !is.na(CRUKID)]
}

plot_nsclc_overview <- function(dt, outfile, genes = NULL) {
  d <- copy(dt)
  if (!is.null(genes) && length(genes) > 0) {
    pat <- paste(genes, collapse = "|")
    d <- d[grepl(pat, Gene)]
  }
  # collapse flags
  d[, flag2 := ifelse(
    grepl("LOH", flag), "LOH",
    ifelse(grepl("LOI", flag), "LOI",
           ifelse(grepl("Imprinted", flag), "Imprinted",
                  ifelse(grepl("no_info", flag), "no_info", "other")))
  )]
  d <- d[flag2 %in% c("Imprinted", "LOI", "LOH")]
  # sample-level unique
  d[, sample_key := paste(CRUKID, sample, sep = "_")]
  sumtab <- d[, .N, by = .(Gene, flag2)]
  sumtab[, Gene := factor(Gene, levels = unique(Gene[order(-N)]))]

  p <- ggplot(sumtab, aes(x = Gene, y = N, fill = flag2)) +
    geom_col(width = 0.75) +
    scale_fill_manual(
      name = NULL,
      values = c(
        Imprinted = "grey70",
        LOI = "firebrick3",
        LOH = "steelblue3"
      )
    ) +
    labs(
      x = NULL,
      y = "Tumour samples",
      title = "S8C - NSCLC imprinted-locus status (TRACERx RRBS)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(size = 10),
      legend.position = "bottom"
    )
  ggsave(outfile, p, width = 7, height = 3.5, units = "in")
  message("Wrote ", outfile)
  invisible(sumtab)
}

plot_sample_strip <- function(dt, gene_pat, outfile, title) {
  d <- copy(dt)[grepl(gene_pat, Gene)]
  d[, flag2 := ifelse(
    grepl("LOH", flag), "LOH",
    ifelse(grepl("LOI", flag), "LOI",
           ifelse(grepl("Imprinted", flag), "Imprinted", "other"))
  )]
  d <- d[flag2 %in% c("Imprinted", "LOI", "LOH")]
  d[, sample_key := paste(CRUKID, sample, sep = "_")]
  # order samples by presence of LOI/LOH
  ord <- d[, .(score = any(flag2 %in% c("LOI", "LOH"))), by = sample_key]
  d[, sample_key := factor(
    sample_key,
    levels = ord[order(-score, sample_key)]$sample_key
  )]
  p <- ggplot(d, aes(x = Gene, y = sample_key, colour = flag2)) +
    geom_point(size = 2) +
    scale_colour_manual(
      values = c(
        Imprinted = "grey55",
        LOI = "firebrick3",
        LOH = "steelblue3"
      )
    ) +
    labs(x = NULL, y = NULL, title = title, colour = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 5),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
      plot.title = element_text(size = 9),
      legend.position = "bottom"
    )
  ggsave(outfile, p, width = 4.5, height = 6, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

################################################################################
# CLI
################################################################################

if (sys.nframe() == 0L) {

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figureS8",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory [default %default]"),
  make_option(c("--tx_dir"), type = "character",
              default = "~/rdeocd/txMeth",
              help = "TRACERx parent (clinical + imprinting RData) [default %default]"),
  make_option(c("--imprint_rdata"), type = "character",
              default = NULL,
              help = "Override imprinting summary RData"),
  make_option(c("--clinical_info"), type = "character",
              default = NULL,
              help = "Clinical TSV for Patient->CRUKID"),
  make_option(c("--panels"), type = "character",
              default = "C",
              help = "Panels to draw (C overview; C_H19 strip) [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)
output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

tx_dir <- path.expand(opt$tx_dir)
clinical_info <- if (!is.null(opt$clinical_info) && nzchar(opt$clinical_info)) {
  path.expand(opt$clinical_info)
} else {
  default_clinical_info(tx_dir)
}
cruk_map <- load_ltx_cruk_map(clinical_info)

rdata <- if (!is.null(opt$imprint_rdata) && nzchar(opt$imprint_rdata)) {
  path.expand(opt$imprint_rdata)
} else {
  file.path(tx_dir, "imprinting_summary_gene_subset_for_plot_Oct2021.RData")
}
if (!file.exists(rdata)) stop("Missing: ", rdata)
message("Loading ", rdata)
e <- new.env(parent = emptyenv())
load(rdata, envir = e)
nm <- ls(e)[1]
dt <- as.data.table(e[[nm]])
if (!"CRUKID" %in% names(dt)) {
  if (!"patient" %in% names(dt)) stop("Need patient or CRUKID column")
  dt[, CRUKID := cruk_map$CRUKID[match(patient, cruk_map$Patient)]]
  dt[is.na(CRUKID), CRUKID := patient]
}
# drop Patient id from published table
out_tab <- copy(dt)
if ("patient" %in% names(out_tab)) out_tab[, patient := NULL]
fwrite(out_tab, file.path(output_dir, "S8C_imprinting_summary_by_sample.tsv"), sep = "\t")

panels <- toupper(trimws(unlist(strsplit(opt$panels, ",", fixed = TRUE))))

if ("C" %in% panels) {
  tab <- plot_nsclc_overview(
    dt,
    file.path(output_dir, "S8C_NSCLC_imprinting_overview.pdf")
  )
  fwrite(tab, file.path(output_dir, "S8C_NSCLC_imprinting_overview_counts.tsv"), sep = "\t")
}

if ("C_H19" %in% panels || "H19" %in% panels) {
  plot_sample_strip(
    dt, "H19",
    file.path(output_dir, "S8C_H19_sample_strip.pdf"),
    title = "S8C - H19 status by sample (CRUKID_Rx)"
  )
}

message(
  "S8A/B locus tracks and S8D CRC overview are not regenerated here; see README."
)
message("Done.")

} # end main
