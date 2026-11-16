#' @title Cohort: within-segment LogR noise pre vs post LogR correction
#' @author Lili Larose Cadieux
#' @description For each tumour sample, use \strong{fixed} ASCAT.m segment
#' boundaries (prefer \code{re_segmented/*.ascat.frag.resegmented.RData};
#' for single-region patients CRUK0010 / CRUK0047 use \code{raw/*.ascat.frag.RData})
#' and score length-weighted mean within-segment SD of raw
#' \code{LogR_Tumour} vs corrected \code{LogR_Tumour_corr} from SNPs tables.
#'
#' BAF is unchanged by LogR correction under fixed segments; this script
#' measures the LogR effect only (no ASPCF re-run).
#'
#' Outputs: per-segment table (top-10 Δ flagged), sample + patient summaries,
#' paired / waterfall plots.

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

resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) {
    return(path)
  }
  file.path(base_dir, path)
}

is_normal_sample <- function(sample_id) {
  grepl("^N([0-9]*|_.*)?$", sample_id)
}

is_tumour_region <- function(sample_id) {
  grepl("^R[0-9]+$", sample_id)
}

#' Single-region patients: no multi-region resegmentation
SINGLE_REGION_PATIENTS <- c("CRUK0010", "CRUK0047")

resolve_frag_path <- function(patient_id, sample_id, sample_dir) {
  reseg <- file.path(
    sample_dir, "re_segmented",
    paste0(patient_id, ".", sample_id, ".ascat.frag.resegmented.RData")
  )
  raw_frag <- file.path(
    sample_dir, "raw",
    paste0(patient_id, ".", sample_id, ".ascat.frag.RData")
  )
  frag <- file.path(
    sample_dir,
    paste0(patient_id, ".", sample_id, ".ascat.frag.RData")
  )
  if (patient_id %in% SINGLE_REGION_PATIENTS) {
    if (file.exists(raw_frag)) {
      return(list(path = raw_frag, source = "raw_frag"))
    }
    if (file.exists(frag)) {
      return(list(path = frag, source = "frag"))
    }
    return(NULL)
  }
  if (file.exists(reseg)) {
    return(list(path = reseg, source = "reseg"))
  }
  # Fallback if reseg missing but raw/frag present
  if (file.exists(raw_frag)) {
    return(list(path = raw_frag, source = "raw_frag_fallback"))
  }
  if (file.exists(frag)) {
    return(list(path = frag, source = "frag_fallback"))
  }
  NULL
}

resolve_snps_path <- function(patient_id, sample_id, sample_dir) {
  txt <- file.path(sample_dir, paste0(patient_id, ".", sample_id, ".SNPs.txt"))
  rdata <- file.path(sample_dir, paste0(patient_id, ".", sample_id, ".SNPs.RData"))
  if (file.exists(txt)) {
    return(list(path = txt, format = "txt"))
  }
  if (file.exists(rdata)) {
    return(list(path = rdata, format = "rdata"))
  }
  NULL
}

load_snps <- function(snps_info) {
  if (snps_info$format == "txt") {
    dt <- fread(snps_info$path)
  } else {
    e <- new.env(parent = emptyenv())
    load(snps_info$path, envir = e)
    dt <- NULL
    for (nm in ls(e)) {
      o <- e[[nm]]
      if (is.data.frame(o) || inherits(o, "data.table")) {
        dt <- as.data.table(o)
        break
      }
    }
    if (is.null(dt)) {
      stop("No data.frame in: ", snps_info$path)
    }
  }
  need <- c("LogR_Tumour", "LogR_Tumour_corr")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("SNPs missing ", paste(miss, collapse = ", "), " in ", snps_info$path)
  }
  dt[]
}

load_frag <- function(frag_path) {
  e <- new.env(parent = emptyenv())
  load(frag_path, envir = e)
  if (!exists("ascat.frag", envir = e)) {
    stop("ascat.frag not in: ", frag_path)
  }
  e$ascat.frag
}

#' Discover tumour regions with SNPs + segment object
discover_samples <- function(tx_dir, patient_filter = NULL) {
  patients <- list.dirs(tx_dir, recursive = FALSE, full.names = FALSE)
  patients <- patients[grepl("^[A-Za-z]+[0-9]", patients)]
  if (!is.null(patient_filter)) {
    patients <- intersect(patients, patient_filter)
  }
  rows <- list()
  for (patient_id in patients) {
    cnd <- file.path(tx_dir, patient_id, "Copy_number")
    if (!dir.exists(cnd)) {
      next
    }
    for (sample_id in list.dirs(cnd, recursive = FALSE, full.names = FALSE)) {
      if (is_normal_sample(sample_id) || !is_tumour_region(sample_id)) {
        next
      }
      sample_dir <- file.path(cnd, sample_id)
      snps <- resolve_snps_path(patient_id, sample_id, sample_dir)
      frag <- resolve_frag_path(patient_id, sample_id, sample_dir)
      if (is.null(snps) || is.null(frag)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.table(
        patient_id = patient_id,
        sample_id = sample_id,
        snps_path = snps$path,
        snps_format = snps$format,
        frag_path = frag$path,
        frag_source = frag$source
      )
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0) {
    stop("No tumour samples with SNPs + frag under: ", tx_dir)
  }
  setorder(out, patient_id, sample_id)
  out[]
}

#' Per-segment within-segment LogR SD (raw vs corr) on fixed fragment segments
segment_logr_noise <- function(logr_raw, logr_corr, logr_seg, chrom,
                               min_n = 10L, autosomes_only = TRUE) {
  stopifnot(length(logr_raw) == length(logr_corr))
  stopifnot(length(logr_raw) == length(logr_seg))
  stopifnot(length(logr_raw) == length(chrom))
  chrom <- as.character(gsub("^chr", "", chrom))
  if (isTRUE(autosomes_only)) {
    keep <- chrom %in% as.character(1:22)
    logr_raw <- logr_raw[keep]
    logr_corr <- logr_corr[keep]
    logr_seg <- logr_seg[keep]
    chrom <- chrom[keep]
  }
  rows <- list()
  seg_id <- 0L
  for (ch in unique(chrom)) {
    ii <- which(chrom == ch)
    r <- rle(as.numeric(logr_seg[ii]))
    pos <- 1L
    for (k in seq_along(r$lengths)) {
      jj <- ii[pos:(pos + r$lengths[k] - 1L)]
      pos <- pos + r$lengths[k]
      seg_id <- seg_id + 1L
      n <- length(jj)
      if (n < min_n) {
        next
      }
      sd_raw <- stats::sd(logr_raw[jj], na.rm = TRUE)
      sd_corr <- stats::sd(logr_corr[jj], na.rm = TRUE)
      if (!is.finite(sd_raw) || !is.finite(sd_corr)) {
        next
      }
      rows[[length(rows) + 1L]] <- data.table(
        seg_id = seg_id,
        chrom = ch,
        n_snp = n,
        seg_logr = r$values[k],
        start_idx = jj[1],
        end_idx = jj[length(jj)],
        noise_raw = sd_raw,
        noise_corr = sd_corr,
        delta = sd_raw - sd_corr
      )
    }
  }
  rbindlist(rows, fill = TRUE)
}

summarise_sample_from_segments <- function(seg_dt) {
  if (nrow(seg_dt) == 0) {
    return(list(
      noise_raw = NA_real_, noise_corr = NA_real_, delta = NA_real_,
      n_seg = 0L, n_snp = 0L
    ))
  }
  w <- seg_dt$n_snp
  list(
    noise_raw = sum(seg_dt$noise_raw * w) / sum(w),
    noise_corr = sum(seg_dt$noise_corr * w) / sum(w),
    delta = sum(seg_dt$delta * w) / sum(w),
    n_seg = nrow(seg_dt),
    n_snp = sum(w)
  )
}

process_sample <- function(row, min_n = 10L, autosomes_only = TRUE) {
  message(
    "=== ", row$patient_id, " ", row$sample_id,
    " [", row$frag_source, "] ==="
  )
  dt <- load_snps(list(path = row$snps_path, format = row$snps_format))
  frag <- load_frag(row$frag_path)
  n <- nrow(frag$Tumor_LogR)
  if (nrow(dt) != n) {
    stop(
      "Row mismatch SNPs (", nrow(dt), ") vs frag (", n, ") for ",
      row$patient_id, " ", row$sample_id
    )
  }
  logr_seg <- as.numeric(frag$Tumor_LogR_segmented[, 1])
  chrom <- as.character(frag$SNPpos[, 1])
  seg_dt <- segment_logr_noise(
    logr_raw = dt$LogR_Tumour,
    logr_corr = dt$LogR_Tumour_corr,
    logr_seg = logr_seg,
    chrom = chrom,
    min_n = min_n,
    autosomes_only = autosomes_only
  )
  seg_dt[, `:=`(
    patient_id = row$patient_id,
    sample_id = row$sample_id,
    frag_source = row$frag_source
  )]
  setcolorder(
    seg_dt,
    c("patient_id", "sample_id", "frag_source", "seg_id", "chrom",
      "n_snp", "seg_logr", "start_idx", "end_idx",
      "noise_raw", "noise_corr", "delta")
  )
  sm <- summarise_sample_from_segments(seg_dt)
  sample_row <- data.table(
    patient_id = row$patient_id,
    sample_id = row$sample_id,
    frag_source = row$frag_source,
    noise_raw = sm$noise_raw,
    noise_corr = sm$noise_corr,
    delta = sm$delta,
    n_seg = sm$n_seg,
    n_snp = sm$n_snp,
    logr_sd_raw_global = stats::sd(dt$LogR_Tumour, na.rm = TRUE),
    logr_sd_corr_global = stats::sd(dt$LogR_Tumour_corr, na.rm = TRUE)
  )
  list(segments = seg_dt, sample = sample_row)
}

summarise_patients <- function(sample_dt) {
  sample_dt[
    is.finite(delta),
    .(
      n_samples = .N,
      noise_raw_mean = mean(noise_raw, na.rm = TRUE),
      noise_corr_mean = mean(noise_corr, na.rm = TRUE),
      delta_mean = mean(delta, na.rm = TRUE),
      delta_se = {
        x <- delta[is.finite(delta)]
        if (length(x) <= 1) NA_real_ else stats::sd(x) / sqrt(length(x))
      },
      frag_sources = paste(unique(frag_source), collapse = ",")
    ),
    by = patient_id
  ]
}

flag_top_deltas <- function(seg_dt, n_top = 10L) {
  d <- copy(seg_dt)
  d[, rank_delta := NA_integer_]
  ord <- order(d$delta, decreasing = TRUE, na.last = TRUE)
  top <- head(ord, n_top)
  d$rank_delta[top] <- seq_along(top)
  d[, top10 := !is.na(rank_delta)]
  d[]
}

#' Patient -> CRUKID map from TRACERx clinical table (as in other figures)
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
  clinical_info <- path.expand(clinical_info)
  if (!file.exists(clinical_info)) {
    stop("Clinical table not found: ", clinical_info)
  }
  dt <- fread(clinical_info)
  need <- c("Patient", "CRUKID")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("Clinical table missing columns: ", paste(miss, collapse = ", "))
  }
  mp <- unique(dt[, .(Patient, CRUKID)])
  mp <- mp[!is.na(Patient) & !is.na(CRUKID) & nzchar(CRUKID)]
  if (anyDuplicated(mp$Patient)) {
    mp <- mp[!duplicated(Patient)]
  }
  mp[]
}

annotate_cruk <- function(dt, ltx_col = "patient_id", map) {
  d <- copy(dt)
  d[, CRUKID := map$CRUKID[match(get(ltx_col), map$Patient)]]
  missing <- unique(d[is.na(CRUKID) | !nzchar(CRUKID), get(ltx_col)])
  if (length(missing) > 0) {
    warning(
      "No CRUKID for: ", paste(missing, collapse = ", "),
      " - falling back to input id"
    )
    d[is.na(CRUKID) | !nzchar(CRUKID), CRUKID := get(ltx_col)]
  }
  d[]
}

#' Publication tables: CRUKID only (drop internal Patient id)
format_output_ids <- function(dt) {
  d <- copy(dt)
  if (!"CRUKID" %in% names(d)) {
    stop("format_output_ids() expects a CRUKID column")
  }
  if ("patient_id" %in% names(d)) {
    d[, patient_id := NULL]
  }
  # Prefer CRUKID as first column
  cols <- c("CRUKID", setdiff(names(d), "CRUKID"))
  setcolorder(d, cols)
  d[]
}

plot_patient_waterfall <- function(patient_dt, outfile, cruk_map = NULL) {
  d <- copy(patient_dt)
  if (!is.null(cruk_map)) {
    d <- annotate_cruk(d, "patient_id", cruk_map)
  } else if (!"CRUKID" %in% names(d)) {
    d[, CRUKID := patient_id]
  }
  d[, CRUKID := factor(CRUKID, levels = CRUKID[order(delta_mean)])]
  cohort_mean <- mean(d$delta_mean, na.rm = TRUE)
  p <- ggplot(d, aes(x = CRUKID, y = delta_mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_hline(
      yintercept = cohort_mean, linewidth = 0.4,
      linetype = "dashed", colour = "firebrick3"
    ) +
    geom_col(width = 0.7, fill = "steelblue3") +
    geom_errorbar(
      aes(
        ymin = delta_mean - ifelse(is.finite(delta_se), delta_se, 0),
        ymax = delta_mean + ifelse(is.finite(delta_se), delta_se, 0)
      ),
      width = 0.15, linewidth = 0.25
    ) +
    labs(
      x = NULL,
      y = expression(Delta ~ "LogR noise (raw - corr)"),
      title = "Patient-mean within-segment LogR SD reduction"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 8, height = 3.5, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

plot_patient_paired <- function(patient_dt, outfile, cruk_map = NULL,
                                n_label = 5L) {
  d0 <- copy(patient_dt)
  if (!is.null(cruk_map)) {
    d0 <- annotate_cruk(d0, "patient_id", cruk_map)
  } else if (!"CRUKID" %in% names(d0)) {
    d0[, CRUKID := patient_id]
  }
  d0 <- d0[is.finite(delta_mean) & is.finite(noise_raw_mean) &
             is.finite(noise_corr_mean)]
  n_label <- min(as.integer(n_label), nrow(d0))
  top_ids <- d0[order(-delta_mean)][seq_len(n_label)]$CRUKID

  d <- melt(
    d0,
    id.vars = c("CRUKID", "delta_mean"),
    measure.vars = c("noise_raw_mean", "noise_corr_mean"),
    variable.name = "which",
    value.name = "noise"
  )
  d[, which := factor(
    which,
    levels = c("noise_raw_mean", "noise_corr_mean"),
    labels = c("Raw LogR", "Corrected LogR")
  )]
  d_top <- d[CRUKID %in% top_ids]
  d_lab <- d_top[which == "Raw LogR"]

  p <- ggplot(d, aes(x = which, y = noise)) +
    geom_boxplot(
      width = 0.45, outlier.shape = NA,
      fill = "grey90", colour = "grey40", linewidth = 0.4
    ) +
    geom_line(
      data = d[!CRUKID %in% top_ids],
      aes(group = CRUKID),
      colour = "grey75", linewidth = 0.25, alpha = 0.7
    ) +
    geom_line(
      data = d_top,
      aes(group = CRUKID),
      colour = "firebrick3", linewidth = 0.5, alpha = 0.9
    ) +
    geom_point(
      data = d[!CRUKID %in% top_ids],
      size = 1.1, alpha = 0.55, colour = "grey40"
    ) +
    geom_point(
      data = d_top,
      size = 1.7, colour = "firebrick3"
    ) +
    geom_text(
      data = d_lab,
      aes(label = CRUKID),
      hjust = 1.08, size = 2.4, colour = "firebrick4",
      check_overlap = TRUE
    ) +
    labs(
      x = NULL,
      y = "Mean within-segment LogR SD",
      title = "LogR noise before / after correction (patient means)",
      subtitle = paste0(
        "Labels: top ", n_label, " largest delta (raw - corr)"
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8, colour = "grey30")
    )
  ggsave(outfile, p, width = 4.5, height = 3.8, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

plot_sample_paired <- function(sample_dt, outfile, n_tie = 10L,
                               cruk_map = NULL) {
  ok <- sample_dt[is.finite(delta) & is.finite(noise_raw) & is.finite(noise_corr)]
  if (!is.null(cruk_map)) {
    ok <- annotate_cruk(ok, "patient_id", cruk_map)
  } else if (!"CRUKID" %in% names(ok)) {
    ok[, CRUKID := patient_id]
  }
  # Publication sample label: CRUKID_Rx (matches figure3FG CRUK_IDs)
  ok[, sample_key := paste(CRUKID, sample_id, sep = "_")]
  n_tie <- min(as.integer(n_tie), nrow(ok))
  top_keys <- ok[order(-delta)][seq_len(n_tie)]$sample_key

  d <- melt(
    ok,
    id.vars = c("patient_id", "sample_id", "CRUKID", "sample_key", "delta"),
    measure.vars = c("noise_raw", "noise_corr"),
    variable.name = "which",
    value.name = "noise"
  )
  d[, which := factor(
    which,
    levels = c("noise_raw", "noise_corr"),
    labels = c("Raw LogR", "Corrected LogR")
  )]
  d_tie <- d[sample_key %in% top_keys]
  # Label top deltas at the raw end
  d_lab <- d_tie[which == "Raw LogR"]

  p <- ggplot(d, aes(x = which, y = noise)) +
    geom_boxplot(
      width = 0.45, outlier.shape = NA,
      fill = "grey90", colour = "grey40", linewidth = 0.4
    ) +
    geom_line(
      data = d_tie,
      aes(group = sample_key),
      colour = "firebrick3", linewidth = 0.45, alpha = 0.85
    ) +
    geom_point(
      data = d[!sample_key %in% top_keys],
      size = 1.0, alpha = 0.35, colour = "grey45",
      position = position_jitter(width = 0.06, height = 0, seed = 1)
    ) +
    geom_point(
      data = d_tie,
      size = 1.6, colour = "firebrick3"
    ) +
    geom_text(
      data = d_lab,
      aes(label = sample_key),
      hjust = 1.05, size = 2.2, colour = "firebrick4",
      check_overlap = TRUE
    ) +
    labs(
      x = NULL,
      y = "Mean within-segment LogR SD",
      title = "LogR noise before / after correction (all samples)",
      subtitle = paste0("Lines: top ", n_tie, " largest delta (raw - corr)")
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8, colour = "grey30")
    )
  ggsave(outfile, p, width = 4.5, height = 3.8, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

################################################################################
# CLI (only when executed as Rscript, not when sourced)
################################################################################

if (sys.nframe() == 0L) {

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
  make_option(c("--patient_id"), type = "character",
              default = NULL,
              help = "Optional comma-separated patient filter"),
  make_option(c("--sample_id"), type = "character",
              default = NULL,
              help = "Optional comma-separated sample filter"),
  make_option(c("--max_samples"), type = "integer",
              default = NULL,
              help = "Optional cap on samples (smoke tests)"),
  make_option(c("--min_snp_n"), type = "integer",
              default = 10L,
              help = "Min SNPs per segment [default %default]"),
  make_option(c("--top_n"), type = "integer",
              default = 10L,
              help = "Flag N largest segment deltas [default %default]"),
  make_option(c("--n_tie"), type = "integer",
              default = 10L,
              help = "N largest sample deltas to connect on paired plot [default %default]"),
  make_option(c("--autosomes_only"), type = "logical",
              default = TRUE,
              help = "Restrict to chr1-22 [default %default]"),
  make_option(c("--clinical_info"), type = "character",
              default = NULL,
              help = paste(
                "Clinical TSV with Patient/CRUKID",
                "(default: tx_dir/tx100clinicalinfo_RRBS_published_status_...)"
              ))
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
clinical_info <- if (!is.null(opt$clinical_info) && nzchar(opt$clinical_info)) {
  path.expand(opt$clinical_info)
} else {
  default_clinical_info(tx_dir)
}
cruk_map <- load_ltx_cruk_map(clinical_info)
message("Loaded Patient->CRUKID map: ", nrow(cruk_map), " patients from ", clinical_info)

patient_filter <- NULL
if (!is.null(opt$patient_id) && nzchar(opt$patient_id)) {
  patient_filter <- trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE)))
}

samples <- discover_samples(tx_dir, patient_filter = patient_filter)
if (!is.null(opt$sample_id) && nzchar(opt$sample_id)) {
  sid <- trimws(unlist(strsplit(opt$sample_id, ",", fixed = TRUE)))
  samples <- samples[sample_id %in% sid]
}
if (!is.null(opt$max_samples) && is.finite(opt$max_samples)) {
  samples <- head(samples, as.integer(opt$max_samples))
}
if (nrow(samples) == 0) {
  stop("No samples left after filters")
}
message("Samples to process: ", nrow(samples))
message(
  "  reseg: ", sum(samples$frag_source == "reseg"),
  "  single-region raw: ",
  sum(samples$frag_source %in% c("raw_frag", "frag")),
  "  fallback: ",
  sum(grepl("fallback", samples$frag_source))
)
print(samples[, .N, by = .(frag_source)])

seg_rows <- vector("list", nrow(samples))
sample_rows <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  res <- tryCatch(
    process_sample(
      samples[i],
      min_n = as.integer(opt$min_snp_n),
      autosomes_only = isTRUE(opt$autosomes_only)
    ),
    error = function(e) {
      warning(
        "Failed ", samples$patient_id[i], " ", samples$sample_id[i],
        ": ", conditionMessage(e)
      )
      list(
        segments = NULL,
        sample = data.table(
          patient_id = samples$patient_id[i],
          sample_id = samples$sample_id[i],
          frag_source = samples$frag_source[i],
          noise_raw = NA_real_,
          noise_corr = NA_real_,
          delta = NA_real_,
          error = conditionMessage(e)
        )
      )
    }
  )
  seg_rows[[i]] <- res$segments
  sample_rows[[i]] <- res$sample
}

seg_dt <- rbindlist(seg_rows, fill = TRUE)
sample_dt <- rbindlist(sample_rows, fill = TRUE)

if (nrow(seg_dt) > 0) {
  seg_dt <- flag_top_deltas(seg_dt, n_top = as.integer(opt$top_n))
}

seg_out <- file.path(output_dir, "LogR_noise_vs_correction_by_segment.tsv")
fwrite(seg_dt, seg_out, sep = "\t")
message("Wrote ", seg_out)

top_out <- file.path(output_dir, "LogR_noise_top10_segment_deltas.tsv")
fwrite(seg_dt[top10 == TRUE][order(rank_delta)], top_out, sep = "\t")
message("Wrote ", top_out)

sample_out <- file.path(output_dir, "LogR_noise_vs_correction_by_sample.tsv")
fwrite(sample_dt, sample_out, sep = "\t")
message("Wrote ", sample_out)

ok <- sample_dt[is.finite(delta)]
if (nrow(ok) == 0) {
  stop("No successful sample results; see ", sample_out)
}

patient_dt <- summarise_patients(ok)
patient_out <- file.path(output_dir, "LogR_noise_vs_correction_by_patient.tsv")
fwrite(patient_dt, patient_out, sep = "\t")
message("Wrote ", patient_out)

# Annotate CRUKID; published tables drop internal Patient id
ok <- annotate_cruk(ok, "patient_id", cruk_map)
patient_dt <- annotate_cruk(patient_dt, "patient_id", cruk_map)
if (nrow(seg_dt) > 0) {
  seg_dt <- annotate_cruk(seg_dt, "patient_id", cruk_map)
  fwrite(format_output_ids(seg_dt), seg_out, sep = "\t")
  fwrite(
    format_output_ids(seg_dt[top10 == TRUE][order(rank_delta)]),
    top_out, sep = "\t"
  )
}
fwrite(format_output_ids(ok), sample_out, sep = "\t")
fwrite(format_output_ids(patient_dt), patient_out, sep = "\t")

plot_sample_paired(
  ok,
  file.path(output_dir, "LogR_noise_paired_by_sample.pdf"),
  n_tie = as.integer(opt$n_tie),
  cruk_map = cruk_map
)
plot_patient_paired(
  patient_dt,
  file.path(output_dir, "LogR_noise_paired_by_patient.pdf"),
  cruk_map = cruk_map,
  n_label = as.integer(opt$n_tie)
)
plot_patient_waterfall(
  patient_dt,
  file.path(output_dir, "LogR_noise_delta_by_patient_waterfall.pdf"),
  cruk_map = cruk_map
)

message(
  "Cohort mean sample delta: ",
  round(mean(ok$delta, na.rm = TRUE), 5),
  " | patient delta: ",
  round(mean(patient_dt$delta_mean, na.rm = TRUE), 5)
)
message("Top-", opt$top_n, " segment deltas:")
print(seg_dt[top10 == TRUE][order(rank_delta),
  .(rank_delta, CRUKID, sample_id, chrom, n_snp, noise_raw, noise_corr, delta)
])
message("Done.")

} # end main (Rscript only)
