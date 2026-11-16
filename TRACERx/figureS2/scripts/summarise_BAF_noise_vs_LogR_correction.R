#' @title Cohort: BAF within-segment noise with vs without LogR correction
#' @author Lili Larose Cadieux
#' @description For each tumour sample with ASCAT.m SNPs + ascat.bc, re-run
#' ASPCF twice (raw \code{LogR_Tumour} vs corrected \code{LogR_Tumour_corr}),
#' measure length-weighted mean within-segment BAF SD among heterozygous SNPs,
#' then average Δ by patient.
#'
#' ASPCF is sourced from a paper-local copy
#' (\code{scripts/lib/ascat.aspcf.modif.R}) patched for R ≥ 4.3
#' (\code{&&} length-1 rules). Pipeline originals are not modified.
#'
#' Outputs:
#' * sample-level TSV (noise_raw, noise_corr, delta, …)
#' * patient-level TSV (mean ± SE)
#' * paired and waterfall PDFs of patient-averaged Δ

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
  grepl("^N[0-9]*$", sample_id)
}

#' Discover tumour samples with SNPs.txt + ascat.bc under tx_dir
discover_samples <- function(tx_dir, patient_filter = NULL) {
  snp_files <- list.files(
    tx_dir, pattern = "\\.SNPs\\.txt$",
    recursive = TRUE, full.names = TRUE
  )
  # Expect .../<patient>/Copy_number/<sample>/<patient>.<sample>.SNPs.txt
  keep <- grepl("/Copy_number/", snp_files)
  snp_files <- snp_files[keep]
  rows <- lapply(snp_files, function(f) {
    bn <- basename(f)
    m <- regmatches(bn, regexec("^([^.]+)\\.([^.]+)\\.SNPs\\.txt$", bn))[[1]]
    if (length(m) < 3) {
      return(NULL)
    }
    patient_id <- m[2]
    sample_id <- m[3]
    if (is_normal_sample(sample_id)) {
      return(NULL)
    }
    if (!is.null(patient_filter) && !(patient_id %in% patient_filter)) {
      return(NULL)
    }
    bc <- file.path(
      dirname(f),
      paste0(patient_id, ".", sample_id, ".ascat.bc.RData")
    )
    if (!file.exists(bc)) {
      return(NULL)
    }
    data.table(
      patient_id = patient_id,
      sample_id = sample_id,
      snps_path = f,
      bc_path = bc
    )
  })
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out) == 0) {
    stop("No tumour samples with SNPs.txt + ascat.bc under: ", tx_dir)
  }
  setorder(out, patient_id, sample_id)
  out[]
}

#' Prepare ascat.bc for ASPCF with a chosen LogR vector
prepare_ascat_for_aspcf <- function(ascat.bc, logr) {
  o <- ascat.bc
  if (length(logr) != nrow(o$Tumor_LogR)) {
    stop(
      "LogR length (", length(logr),
      ") != Tumor_LogR rows (", nrow(o$Tumor_LogR), ")"
    )
  }
  # Keep structure; overwrite tumour LogR
  if (is.data.frame(o$Tumor_LogR)) {
    o$Tumor_LogR[[1]] <- as.numeric(logr)
  } else {
    o$Tumor_LogR[, 1] <- as.numeric(logr)
  }
  o$Tumor_LogR_segmented <- NULL
  o$Tumor_BAF_segmented <- NULL
  # ascat.aspcf.modif expects $chrom (= chromosome-arm index lists)
  o$chrom <- o$chr
  if (is.null(o$sexchromosomes) && !is.null(o$sex.chromosomes)) {
    o$sexchromosomes <- o$sex.chromosomes
  }
  if (is.null(colnames(o$Tumor_LogR))) {
    colnames(o$Tumor_LogR) <- o$samples
  }
  if (is.null(colnames(o$Tumor_BAF))) {
    colnames(o$Tumor_BAF) <- o$samples
  }
  o
}

#' Length-weighted mean within-segment BAF SD (heterozygous probes only)
within_segment_baf_noise <- function(ascat_frag, min_n = 10L, sample_idx = 1L) {
  bs <- ascat_frag$Tumor_BAF_segmented[[sample_idx]]
  if (is.null(bs) || length(bs) == 0) {
    return(list(
      noise = NA_real_, n_het = 0L, n_seg = 0L, n_seg_used = 0L
    ))
  }
  idx <- as.integer(rownames(bs))
  if (anyNA(idx) || any(idx < 1L) || any(idx > nrow(ascat_frag$Tumor_BAF))) {
    # Fall back to name match if rownames are probe IDs
    idx <- match(rownames(bs), rownames(ascat_frag$Tumor_BAF))
  }
  ok <- !is.na(idx)
  idx <- idx[ok]
  baf <- as.numeric(ascat_frag$Tumor_BAF[idx, sample_idx])
  seg <- as.numeric(bs[ok, 1])
  chrom <- as.character(ascat_frag$SNPpos[idx, 1])

  noise_w <- 0
  wtot <- 0
  n_seg_used <- 0L
  n_seg <- 0L
  for (ch in unique(chrom)) {
    ii <- which(chrom == ch)
    r <- rle(seg[ii])
    pos <- 1L
    for (k in seq_along(r$lengths)) {
      jj <- ii[pos:(pos + r$lengths[k] - 1L)]
      pos <- pos + r$lengths[k]
      n_seg <- n_seg + 1L
      if (length(jj) < min_n) {
        next
      }
      sdv <- stats::sd(baf[jj], na.rm = TRUE)
      if (!is.finite(sdv)) {
        next
      }
      w <- length(jj)
      noise_w <- noise_w + sdv * w
      wtot <- wtot + w
      n_seg_used <- n_seg_used + 1L
    }
  }
  list(
    noise = if (wtot > 0) noise_w / wtot else NA_real_,
    n_het = length(baf),
    n_seg = n_seg,
    n_seg_used = n_seg_used
  )
}

#' Run ASPCF in a temp directory (writes *.PCFed.txt sidecars)
run_aspcf <- function(ascat.bc, logr, penalty = 200) {
  o <- prepare_ascat_for_aspcf(ascat.bc, logr)
  td <- tempfile("aspcf_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  owd <- setwd(td)
  on.exit(setwd(owd), add = TRUE)
  ascat.aspcf.modif(o, selectsamples = 1, penalty = penalty)
}

process_sample <- function(patient_id, sample_id, snps_path, bc_path,
                           penalty = 200, min_n = 10L) {
  message("=== ", patient_id, " ", sample_id, " ===")
  dt <- fread(snps_path)
  need <- c("LogR_Tumour", "LogR_Tumour_corr")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("SNPs missing columns: ", paste(miss, collapse = ", "), " in ", snps_path)
  }
  e <- new.env(parent = emptyenv())
  load(bc_path, envir = e)
  if (!exists("ascat.bc", envir = e)) {
    stop("ascat.bc not found in: ", bc_path)
  }
  ascat.bc <- e$ascat.bc
  if (nrow(dt) != nrow(ascat.bc$Tumor_LogR)) {
    stop(
      "Row mismatch SNPs (", nrow(dt), ") vs ascat.bc (",
      nrow(ascat.bc$Tumor_LogR), ") for ", patient_id, " ", sample_id
    )
  }

  logr_sd_raw <- stats::sd(dt$LogR_Tumour, na.rm = TRUE)
  logr_sd_corr <- stats::sd(dt$LogR_Tumour_corr, na.rm = TRUE)

  message("  ASPCF raw LogR…")
  t0 <- Sys.time()
  frag_raw <- run_aspcf(ascat.bc, dt$LogR_Tumour, penalty = penalty)
  message("    elapsed ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s")
  n_raw <- within_segment_baf_noise(frag_raw, min_n = min_n)

  message("  ASPCF corrected LogR…")
  t1 <- Sys.time()
  frag_corr <- run_aspcf(ascat.bc, dt$LogR_Tumour_corr, penalty = penalty)
  message("    elapsed ", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), "s")
  n_corr <- within_segment_baf_noise(frag_corr, min_n = min_n)

  data.table(
    patient_id = patient_id,
    sample_id = sample_id,
    noise_raw = n_raw$noise,
    noise_corr = n_corr$noise,
    delta = n_raw$noise - n_corr$noise,
    n_het_raw = n_raw$n_het,
    n_het_corr = n_corr$n_het,
    n_seg_raw = n_raw$n_seg,
    n_seg_corr = n_corr$n_seg,
    n_seg_used_raw = n_raw$n_seg_used,
    n_seg_used_corr = n_corr$n_seg_used,
    logr_sd_raw = logr_sd_raw,
    logr_sd_corr = logr_sd_corr,
    logr_sd_delta = logr_sd_raw - logr_sd_corr
  )
}

summarise_patients <- function(sample_dt) {
  sample_dt[
    ,
    .(
      n_samples = .N,
      noise_raw_mean = mean(noise_raw, na.rm = TRUE),
      noise_corr_mean = mean(noise_corr, na.rm = TRUE),
      delta_mean = mean(delta, na.rm = TRUE),
      delta_se = {
        x <- delta[is.finite(delta)]
        if (length(x) <= 1) NA_real_ else stats::sd(x) / sqrt(length(x))
      },
      logr_sd_delta_mean = mean(logr_sd_delta, na.rm = TRUE)
    ),
    by = patient_id
  ]
}

plot_patient_delta_waterfall <- function(patient_dt, outfile) {
  d <- copy(patient_dt)
  d[, patient_id := factor(patient_id, levels = patient_id[order(delta_mean)])]
  cohort_mean <- mean(d$delta_mean, na.rm = TRUE)
  p <- ggplot(d, aes(x = patient_id, y = delta_mean)) +
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
      width = 0.2, linewidth = 0.3
    ) +
    labs(
      x = NULL,
      y = expression(Delta ~ "BAF noise (raw - corr)"),
      title = "Patient-mean within-segment BAF SD reduction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 4.5, height = 3.2, units = "in")
  message("Wrote ", outfile)
  invisible(p)
}

plot_patient_paired <- function(patient_dt, outfile) {
  d <- melt(
    patient_dt,
    id.vars = "patient_id",
    measure.vars = c("noise_raw_mean", "noise_corr_mean"),
    variable.name = "which",
    value.name = "noise"
  )
  d[, which := factor(
    which,
    levels = c("noise_raw_mean", "noise_corr_mean"),
    labels = c("Raw LogR", "Corrected LogR")
  )]
  p <- ggplot(d, aes(x = which, y = noise, group = patient_id)) +
    geom_line(colour = "grey55", linewidth = 0.4) +
    geom_point(aes(colour = patient_id), size = 2) +
    labs(
      x = NULL,
      y = "Mean within-segment BAF SD",
      colour = "Patient",
      title = "BAF noise before / after LogR correction"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10)
    )
  ggsave(outfile, p, width = 4.8, height = 3.2, units = "in")
  message("Wrote ", outfile)
  invisible(p)
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
  make_option(c("--aspcf_script"), type = "character",
              default = "scripts/lib/ascat.aspcf.modif.R",
              help = "Paper-local ASPCF source [default %default]"),
  make_option(c("--patient_id"), type = "character",
              default = NULL,
              help = "Optional comma-separated patient filter (default: all)"),
  make_option(c("--sample_id"), type = "character",
              default = NULL,
              help = "Optional comma-separated sample filter (default: all)"),
  make_option(c("--max_samples"), type = "integer",
              default = NULL,
              help = "Optional cap on number of samples (smoke tests)"),
  make_option(c("--penalty"), type = "integer",
              default = 200L,
              help = "ASPCF penalty [default %default]"),
  make_option(c("--min_het_n"), type = "integer",
              default = 10L,
              help = "Min heterozygous SNPs per segment [default %default]")
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

aspcf_script <- resolve_path(opt$aspcf_script, work_dir)
if (!file.exists(aspcf_script)) {
  stop("ASPCF script not found: ", aspcf_script)
}
source(aspcf_script)
if (!exists("ascat.aspcf.modif")) {
  stop("ascat.aspcf.modif() not defined after sourcing: ", aspcf_script)
}

tx_dir <- path.expand(opt$tx_dir)
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
print(samples[, .(patient_id, sample_id)])

rows <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  rows[[i]] <- tryCatch(
    process_sample(
      patient_id = samples$patient_id[i],
      sample_id = samples$sample_id[i],
      snps_path = samples$snps_path[i],
      bc_path = samples$bc_path[i],
      penalty = as.integer(opt$penalty),
      min_n = as.integer(opt$min_het_n)
    ),
    error = function(e) {
      warning(
        "Failed ", samples$patient_id[i], " ", samples$sample_id[i],
        ": ", conditionMessage(e)
      )
      data.table(
        patient_id = samples$patient_id[i],
        sample_id = samples$sample_id[i],
        noise_raw = NA_real_,
        noise_corr = NA_real_,
        delta = NA_real_,
        error = conditionMessage(e)
      )
    }
  )
}
sample_dt <- rbindlist(rows, fill = TRUE)

sample_out <- file.path(output_dir, "BAF_noise_vs_LogR_correction_by_sample.tsv")
fwrite(sample_dt, sample_out, sep = "\t")
message("Wrote ", sample_out)

ok <- sample_dt[is.finite(delta)]
if (nrow(ok) == 0) {
  stop("No successful sample results; see ", sample_out)
}

patient_dt <- summarise_patients(ok)
patient_out <- file.path(output_dir, "BAF_noise_vs_LogR_correction_by_patient.tsv")
fwrite(patient_dt, patient_out, sep = "\t")
message("Wrote ", patient_out)

plot_patient_delta_waterfall(
  patient_dt,
  file.path(output_dir, "BAF_noise_delta_by_patient_waterfall.pdf")
)
plot_patient_paired(
  patient_dt,
  file.path(output_dir, "BAF_noise_paired_by_patient.pdf")
)

message(
  "Cohort mean Δ (patient-level): ",
  round(mean(patient_dt$delta_mean, na.rm = TRUE), 5)
)
message("Done.")
