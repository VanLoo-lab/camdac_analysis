#' @title Supplementary Figure S2A–D: RRBS LogR biases and correction
#' @author Lili Larose Cadieux
#' @description Visualise MspI fragment length, GC content and replication
#' timing biases on tumour LogR, and the combined natural-spline correction
#' used by ASCAT.m (\code{LogR_correction} in \code{run_ASCAT_RRBS_cluster_v3.R}).
#'
#' Panels:
#' * A — predicted LogR bias vs log10(MspI fragment length)
#' * B — predicted LogR bias vs fragment GC content
#' * C — predicted LogR bias vs replication timing
#' * D — combined correction: raw vs corrected LogR (+ fitted bias)
#'
#' S2E–G (BAF rule schematics) are not generated here.
#'
#' Extracted from \code{archive/Notes_correcting_logR_Aug2018_Mar2019.R}
#' (NOTES update 02/04/2019) and \code{LogR_correction()}.

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
x <- c("data.table", "ggplot2", "splines", "scales", "optparse", "fst")
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

resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) {
    return(path)
  }
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

#' Resolve CLI patient id (CRUKID (or Patient)) -> list(ltx, cruk)
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

save_panel <- function(plot, outfile, width = 3, height = 3) {
  ggsave(
    filename = outfile, plot = plot, device = "pdf",
    width = width, height = height, units = "in"
  )
  message("Wrote ", outfile)
  invisible(plot)
}

#' Load ASCAT.m SNP table written after LogR_correction
load_snps <- function(tx_dir, patient_id, sample_id) {
  path_txt <- file.path(
    tx_dir, patient_id, "Copy_number", sample_id,
    paste0(patient_id, ".", sample_id, ".SNPs.txt")
  )
  path_rdata <- file.path(
    tx_dir, patient_id, "Copy_number", sample_id,
    paste0(patient_id, ".", sample_id, ".SNPs.RData")
  )
  if (file.exists(path_txt)) {
    dt <- fread(path_txt)
  } else if (file.exists(path_rdata)) {
    e <- new.env(parent = emptyenv())
    load(path_rdata, envir = e)
    # Prefer a data.frame / data.table object
    nms <- ls(e)
    dt <- NULL
    for (nm in nms) {
      o <- e[[nm]]
      if (is.data.frame(o) || inherits(o, "data.table")) {
        dt <- as.data.table(o)
        break
      }
    }
    if (is.null(dt)) {
      stop("No data.frame/data.table in: ", path_rdata)
    }
  } else {
    stop("Missing SNPs table:\n  ", path_txt, "\n  ", path_rdata)
  }
  req <- c("LogR_Tumour", "Log_msp1_length", "GC_content")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) {
    stop("SNPs table missing columns: ", paste(miss, collapse = ", "))
  }
  if (!"LogR_Tumour_corr" %in% names(dt)) {
    warning("LogR_Tumour_corr missing; panel D will refit the joint model")
  }
  if (!"chrom" %in% names(dt) && "CHR" %in% names(dt)) {
    dt[, chrom := gsub("^chr", "", CHR)]
  }
  if ("chrom" %in% names(dt)) {
    dt[, chrom := as.character(gsub("^chr", "", chrom))]
  }
  dt[]
}

#' Annotate SNPs with the replication-timing track best correlated with LogR
#' (mirrors \code{LogR_correction} max |cor| track selection).
annotate_replic <- function(dt, replic_prefix, n_cores = 1L) {
  if ("replic" %in% names(dt) && dt[, any(is.finite(replic))]) {
    return(dt)
  }
  if (!all(c("chrom", "POS", "LogR_Tumour") %in% names(dt))) {
    stop("Need chrom, POS, LogR_Tumour to annotate replication timing")
  }
  chroms <- intersect(as.character(1:22), unique(dt$chrom))
  if (length(chroms) == 0) {
    stop("No autosomal chrom labels found for replication annotation")
  }

  load_chr <- function(chr) {
    f <- paste0(replic_prefix, chr, ".fst")
    if (!file.exists(f)) {
      return(NULL)
    }
    x <- as.data.table(read_fst(f))
    # Expected columns: chr, pos, then multiple cell-type tracks
    pos_col <- if ("pos" %in% names(x)) "pos" else names(x)[2]
    setnames(x, pos_col, "pos", skip_absent = TRUE)
    x[, chrom := as.character(chr)]
    x[]
  }

  if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    replic_list <- parallel::mclapply(chroms, load_chr, mc.cores = n_cores)
  } else {
    replic_list <- lapply(chroms, load_chr)
  }
  replic_list <- Filter(Negate(is.null), replic_list)
  if (length(replic_list) == 0) {
    stop("No replication-timing .fst files found under prefix: ", replic_prefix)
  }
  replic <- rbindlist(replic_list, fill = TRUE)
  track_cols <- setdiff(names(replic), c("chr", "chrom", "pos"))
  track_cols <- track_cols[vapply(
    track_cols, function(cn) is.numeric(replic[[cn]]), logical(1)
  )]
  if (length(track_cols) == 0) {
    stop("No numeric replication tracks in fst files")
  }

  # Nearest replic locus per SNP (rolling join on pos within chrom)
  snp <- dt[, .(chrom = as.character(chrom), POS = as.numeric(POS),
                LogR_Tumour, row_id = .I)]
  setkey(snp, chrom, POS)
  replic[, pos := as.numeric(pos)]
  setkey(replic, chrom, pos)
  joined <- replic[snp, roll = "nearest"]

  cors <- vapply(track_cols, function(cn) {
    suppressWarnings(cor(joined$LogR_Tumour, joined[[cn]], use = "complete.obs"))
  }, numeric(1))
  best <- names(which.max(abs(cors)))
  message("Selected replication track: ", best,
          " (|cor|=", signif(abs(cors[best]), 3), ")")
  dt <- copy(dt)
  dt[, replic := joined[[best]][order(joined$row_id)]]
  dt[]
}

#' Univariate natural-spline bias curve (Notes S2A–C style)
plot_univariate_bias <- function(y, x, xlab, colour, outfile,
                                 xlim = NULL, ylim = c(-1.2, 1.2),
                                 log10_x = FALSE, df_spline = 5L,
                                 n_plot = 1000L, width = 3, height = 3) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]
  x <- x[ok]
  model <- lm(
    y ~ splines::ns(x = x, df = df_spline, intercept = TRUE),
    data = data.frame(y = y, x = x)
  )
  ord <- order(x)
  # Subsample along ordered x for a smooth thin line (match Notes)
  idx <- unique(as.integer(
    seq(1, length(y), length.out = min(n_plot, length(y)))
  ))
  pred <- predict(model)[ord][idx]
  x_plot <- x[ord][idx]
  if (is.null(xlim)) {
    xlim <- range(x_plot, na.rm = TRUE)
  }
  p <- ggplot(data.frame(x = x_plot, y = pred)) +
    geom_line(aes(x = x, y = y), colour = colour, linewidth = 0.8) +
    theme_classic() +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    ylab("Predicted tumour LogR bias") +
    xlab(xlab)
  if (isTRUE(log10_x)) {
    p <- p + scale_x_continuous(labels = math_format(10^.x))
  }
  save_panel(p, outfile, width = width, height = height)
  invisible(list(model = model, plot = p))
}

#' Panel D: combined correction summary
plot_combined_correction <- function(dt, outfile, df_spline = 5L,
                                     width = 3.5, height = 3) {
  dt <- dt[is.finite(LogR_Tumour) & is.finite(Log_msp1_length) &
             is.finite(GC_content) & is.finite(replic)]
  y <- dt$LogR_Tumour
  model <- lm(
    y ~ splines::ns(x = Log_msp1_length, df = df_spline, intercept = TRUE) +
      splines::ns(x = GC_content, df = df_spline, intercept = TRUE) +
      splines::ns(x = replic, df = df_spline, intercept = TRUE),
    data = dt
  )
  dt[, fitted_bias := as.numeric(predict(model))]
  if ("LogR_Tumour_corr" %in% names(dt)) {
    dt[, LogR_corr_plot := LogR_Tumour_corr]
  } else {
    dt[, LogR_corr_plot := LogR_Tumour - fitted_bias]
  }

  long <- rbind(
    dt[, .(LogR = LogR_Tumour, status = "Raw")],
    dt[, .(LogR = LogR_corr_plot, status = "Corrected")]
  )
  long[, status := factor(status, levels = c("Raw", "Corrected"))]

  p <- ggplot(long, aes(x = LogR, colour = status)) +
    theme_classic() +
    geom_density(adjust = 1.2, linewidth = 0.7) +
    scale_colour_manual(
      values = c(Raw = "grey50", Corrected = "firebrick3"),
      name = NULL
    ) +
    coord_cartesian(xlim = c(-3, 3)) +
    xlab("Tumour LogR") +
    ylab("Density") +
    theme(legend.position = "top")
  save_panel(p, outfile, width = width, height = height)

  # Companion: variance summary sidecar (same stem)
  stats <- data.table(
    metric = c("sd_raw", "sd_corrected", "var_explained_R2"),
    value = c(
      sd(dt$LogR_Tumour, na.rm = TRUE),
      sd(dt$LogR_corr_plot, na.rm = TRUE),
      summary(model)$r.squared
    )
  )
  fwrite(stats, sub("\\.pdf$", "_stats.tsv", outfile), sep = "\t")
  message("Wrote ", sub("\\.pdf$", "_stats.tsv", outfile))
  invisible(list(model = model, plot = p, stats = stats))
}

################################################################################
# CLI
################################################################################

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figureS2",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory [default %default]"),
  make_option(c("--tx_dir"), type = "character",
              default = "~/rdeocd/txMeth",
              help = "TRACERx methylation parent directory [default %default]"),
  make_option(c("--pipeline_dir"), type = "character",
              default = "~/rdeocd/pipeline_files",
              help = "CAMDAC pipeline_files directory [default %default]"),
  make_option(c("--patient_id"), type = "character",
              default = "CRUK0062",
              help = "Patient ID as CRUKID (or Patient) (outputs use CRUKID) [default %default]"),
  make_option(c("--clinical_info"), type = "character",
              default = NULL,
              help = "Optional clinical TSV for Patient<->CRUKID map"),
  make_option(c("--sample_id"), type = "character",
              default = "R3",
              help = "Sample / region ID [default %default]"),
  make_option(c("--panels"), type = "character",
              default = "A,B,C,D",
              help = "Comma-separated panels [default %default]"),
  make_option(c("--df_spline"), type = "integer",
              default = 5L,
              help = "Natural spline df (match LogR_correction) [default %default]"),
  make_option(c("--n_cores"), type = "integer",
              default = 4L,
              help = "Cores for reading replic fst files [default %default]")
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
pipeline_dir <- path.expand(opt$pipeline_dir)
ids <- resolve_patient_ids(opt$patient_id, tx_dir, opt$clinical_info)
patient_id <- ids$ltx
cruk_id <- ids$cruk
sample_id <- opt$sample_id
panels <- toupper(trimws(unlist(strsplit(opt$panels, ",", fixed = TRUE))))
df_spline <- as.integer(opt$df_spline)
prefix <- paste0(cruk_id, "_", sample_id)

message("Loading SNPs: ", patient_id, " ", sample_id, " (output ", cruk_id, ")")
dt <- load_snps(tx_dir, patient_id, sample_id)

if (any(panels %in% c("C", "D"))) {
  replic_prefix <- file.path(
    pipeline_dir, "1000genomes_2012_v3_repliTiming",
    "1000_genomes_replication_timing_chr_"
  )
  message("Annotating replication timing…")
  dt <- annotate_replic(dt, replic_prefix, n_cores = opt$n_cores)
}

################################################################################
# Panels
################################################################################

if ("A" %in% panels) {
  q <- quantile(dt$Log_msp1_length, c(0.05, 0.95), na.rm = TRUE)
  plot_univariate_bias(
    y = dt$LogR_Tumour, x = dt$Log_msp1_length,
    xlab = "Log Msp1 fragment length",
    colour = "cornflowerblue",
    outfile = file.path(output_dir, paste0(prefix, "_S2A_Fragment_Length_LogR_bias.pdf")),
    xlim = as.numeric(q), ylim = c(-1.2, 1.2),
    log10_x = TRUE, df_spline = df_spline
  )
}

if ("B" %in% panels) {
  plot_univariate_bias(
    y = dt$LogR_Tumour, x = dt$GC_content,
    xlab = "Fragment GC content",
    colour = "mediumpurple3",
    outfile = file.path(output_dir, paste0(prefix, "_S2B_GC_LogR_bias.pdf")),
    xlim = c(0, 1), ylim = c(-1, 1),
    log10_x = FALSE, df_spline = df_spline
  )
}

if ("C" %in% panels) {
  plot_univariate_bias(
    y = dt$LogR_Tumour, x = dt$replic,
    xlab = "Replication timing",
    colour = "orange3",
    outfile = file.path(output_dir, paste0(prefix, "_S2C_replic_LogR_bias.pdf")),
    xlim = NULL, ylim = c(-0.3, 0),
    log10_x = FALSE, df_spline = df_spline
  )
}

if ("D" %in% panels) {
  plot_combined_correction(
    dt,
    outfile = file.path(output_dir, paste0(prefix, "_S2D_LogR_raw_vs_corrected.pdf")),
    df_spline = df_spline
  )
}

message("Done.")
