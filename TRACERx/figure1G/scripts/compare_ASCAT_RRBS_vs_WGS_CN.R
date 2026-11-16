#' @title Figure 1G: RRBS ASCAT.m vs WGS allele-specific copy number
#' @author Lili Larose Cadieux
#' @description Direct comparison of ASCAT(.m) and WGS (ASCAT / Battenberg-style)
#' allele-specific copy number estimates derived from RRBS and WGS.
#'
#' Output is a segmented copy-number profile rounded to integer copy number.
#' Colouring (total = nA+nB, minor = nB; matches the paper / Notes CN block):
#' * agree — minor allele pale grey, total copy number dark grey
#' * disagree (minor) — yellow (RRBS), baby blue (WGS)
#' * disagree (total) — orange (RRBS), dark blue (WGS)
#'
#' Extracted from
#' \code{Notes_RRBS_WGS_Genotype_LogR_BAF_comparison.R} (CN comparison block).

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
x <- c("data.table", "ggplot2", "GenomicRanges", "optparse")
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

#' Cap copy-number values for plotting (match legacy Notes behaviour)
cap_cn <- function(x) {
  ifelse(x > 8, 8.3, ifelse(x == 0, 0.2, x))
}

#' Normalise chromosome labels to 1..22, X->23 for ordering / offsets
chrom_to_num <- function(chrom) {
  chrom <- as.character(chrom)
  chrom[chrom == "X"] <- "23"
  chrom[chrom == "chrX"] <- "23"
  chrom <- gsub("^chr", "", chrom)
  as.integer(chrom)
}

#' Load ASCAT(.m) / Battenberg-style segments and return nA/nB data.table
#' @param rdata_path Path to ascat.output*.RData
#' @param platform Label for provenance ("RRBS" or "WGS")
load_ascat_segments <- function(rdata_path, platform) {
  if (!file.exists(rdata_path)) {
    stop("Missing ASCAT output: ", rdata_path)
  }
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  if (!exists("ascat.output", envir = e)) {
    stop("Object 'ascat.output' not found in: ", rdata_path)
  }
  ascat.output <- e$ascat.output

  seg <- as.data.frame(ascat.output$segments, stringsAsFactors = FALSE)
  # ASCAT uses nA/nB; Battenberg-style uses nMajor/nMinor
  if ("nMajor" %in% colnames(seg)) {
    colnames(seg)[colnames(seg) == "nMajor"] <- "nA"
  }
  if ("nMinor" %in% colnames(seg)) {
    colnames(seg)[colnames(seg) == "nMinor"] <- "nB"
  }
  if (!all(c("nA", "nB") %in% colnames(seg))) {
    # fall back to columns 7:8 as in the Notes script
    if (ncol(seg) < 8) {
      stop("Cannot find nA/nB (or nMajor/nMinor) in segments for ", platform)
    }
    colnames(seg)[7:8] <- c("nA", "nB")
  }

  chr_col <- if ("chr" %in% colnames(seg)) {
    "chr"
  } else if ("chromosome" %in% colnames(seg)) {
    "chromosome"
  } else {
    colnames(seg)[1]
  }
  start_col <- if ("startpos" %in% colnames(seg)) "startpos" else "start"
  end_col <- if ("endpos" %in% colnames(seg)) "endpos" else "end"

  dt <- data.table(
    chrom = as.character(seg[[chr_col]]),
    start = as.numeric(seg[[start_col]]),
    end = as.numeric(seg[[end_col]]),
    nA = as.numeric(round(seg$nA)),
    nB = as.numeric(round(seg$nB)),
    platform = platform
  )
  dt[, chrom_num := chrom_to_num(chrom)]
  dt <- dt[!is.na(chrom_num) & chrom_num >= 1 & chrom_num <= 23]
  dt[]
}

#' Build chromosome cumulative offsets for genome-wide x-axis
chrom_offsets <- function(dt) {
  chr <- dt[, .(chr_end = max(end)), by = .(chrom_num)][order(chrom_num)]
  chr[, offset := cumsum(c(0, chr_end[-.N]))]
  chr[, chr_delim := offset + chr_end]
  chr[, chr_label := offset + chr_end / 2]
  chr[]
}

#' Overlap RRBS and WGS segments and assign agreement colours
#' @return data.table ready for ggplot geom_rect
build_agreement_plot_dt <- function(dt_rrbs, dt_wgs) {
  gr_r <- makeGRangesFromDataFrame(
    dt_rrbs, seqnames.field = "chrom_num",
    start.field = "start", end.field = "end",
    keep.extra.columns = TRUE
  )
  gr_w <- makeGRangesFromDataFrame(
    dt_wgs, seqnames.field = "chrom_num",
    start.field = "start", end.field = "end",
    keep.extra.columns = TRUE
  )

  ov <- findOverlaps(gr_r, gr_w)
  if (length(ov) == 0) {
    stop("No overlapping RRBS / WGS copy-number segments")
  }

  inter <- pintersect(gr_r[queryHits(ov)], gr_w[subjectHits(ov)])
  dt_ov <- data.table(
    chrom_num = as.integer(as.character(seqnames(inter))),
    start = start(inter),
    end = end(inter),
    nA_rrbs = gr_r$nA[queryHits(ov)],
    nB_rrbs = gr_r$nB[queryHits(ov)],
    nA_wgs = gr_w$nA[subjectHits(ov)],
    nB_wgs = gr_w$nB[subjectHits(ov)]
  )
  dt_ov <- dt_ov[end > start]

  dt_ov[, tot_rrbs := nA_rrbs + nB_rrbs]
  dt_ov[, tot_wgs := nA_wgs + nB_wgs]
  dt_ov[, agree_total := tot_rrbs == tot_wgs]
  dt_ov[, agree_minor := nB_rrbs == nB_wgs]

  # Plot values (capped), preserving integer identity for agreement
  dt_ov[, `:=`(
    tot_rrbs_plot = cap_cn(tot_rrbs),
    nB_rrbs_plot = cap_cn(nB_rrbs),
    tot_wgs_plot = cap_cn(tot_wgs),
    nB_wgs_plot = cap_cn(nB_wgs)
  )]

  # Long format: one row per series x platform band.
  # WGS slightly below integer, RRBS slightly above, so both colours show on
  # disagreement; identical greys stack when they agree.
  total <- rbind(
    dt_ov[, .(
      chrom_num, start, end,
      ymin = tot_wgs_plot - 0.18, ymax = tot_wgs_plot,
      fill = ifelse(agree_total, "total_agree", "total_wgs_disagree")
    )],
    dt_ov[, .(
      chrom_num, start, end,
      ymin = tot_rrbs_plot, ymax = tot_rrbs_plot + 0.18,
      fill = ifelse(agree_total, "total_agree", "total_rrbs_disagree")
    )]
  )
  minor <- rbind(
    dt_ov[, .(
      chrom_num, start, end,
      ymin = nB_wgs_plot - 0.18, ymax = nB_wgs_plot,
      fill = ifelse(agree_minor, "minor_agree", "minor_wgs_disagree")
    )],
    dt_ov[, .(
      chrom_num, start, end,
      ymin = nB_rrbs_plot, ymax = nB_rrbs_plot + 0.18,
      fill = ifelse(agree_minor, "minor_agree", "minor_rrbs_disagree")
    )]
  )
  rbind(total, minor)
}

#' Plot rounded ASCN comparison for one patient / region
plot_ascat_vs_wgs_cn <- function(dt_plot, chr, outfile,
                                 width = 9.5, height = 1.75) {
  fill_values <- c(
    "minor_agree" = "grey80",           # pale grey
    "total_agree" = "grey30",           # dark grey
    "minor_rrbs_disagree" = "yellow",
    "minor_wgs_disagree" = "lightblue", # baby blue
    "total_rrbs_disagree" = "orange",
    "total_wgs_disagree" = "darkblue"
  )
  fill_labels <- c(
    "minor_agree" = "Minor (both)",
    "total_agree" = "Total (both)",
    "minor_rrbs_disagree" = "Minor (RRBS only)",
    "minor_wgs_disagree" = "Minor (WGS only)",
    "total_rrbs_disagree" = "Total (RRBS only)",
    "total_wgs_disagree" = "Total (WGS only)"
  )

  p <- ggplot() +
    theme_classic() +
    scale_y_continuous(
      name = "Copy number",
      limits = c(0, 9),
      breaks = c(seq(0, 8, 1), 8.5),
      labels = c(seq(0, 8, 1), ">8")
    ) +
    geom_rect(
      data = dt_plot,
      aes(
        xmin = start + offset, xmax = end + offset,
        ymin = ymin, ymax = ymax, fill = fill
      ),
      alpha = 0.85, colour = NA
    ) +
    scale_fill_manual(
      name = "Allele-specific\ncopy number",
      values = fill_values,
      labels = fill_labels,
      breaks = names(fill_values)
    ) +
    geom_vline(
      data = chr, aes(xintercept = chr_delim),
      colour = "darkgrey", alpha = 0.5
    ) +
    geom_vline(xintercept = 0, colour = "darkgrey", alpha = 0.5) +
    geom_text(
      data = chr,
      aes(
        x = chr_label, y = 8.7,
        label = ifelse(chrom_num == 23, "X", chrom_num)
      ),
      size = 3
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  ggsave(
    filename = outfile, plot = p, device = "pdf",
    width = width, height = height, units = "in"
  )
  invisible(p)
}

################################################################################
# Parse command line arguments
################################################################################

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figure1G",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory relative to --work_dir or absolute [default %default]"),
  make_option(c("--tx_dir"), type = "character",
              default = "~/rdeocd/txMeth",
              help = "TRACERx methylation parent directory [default %default]"),
  make_option(c("--patient_id"), type = "character",
              default = "CRUK0062",
              help = "Comma-separated CRUKID(s) (or Patient); outputs use CRUKID [default %default]"),
  make_option(c("--clinical_info"), type = "character",
              default = NULL,
              help = "Optional clinical TSV for Patient<->CRUKID map"),
  make_option(c("--sample_id"), type = "character",
              default = NULL,
              help = paste0(
                "Comma-separated sample / region ID(s), applied to each patient. ",
                "Default: all tumour regions with both RRBS re-segmented and WGS ASCAT outputs"
              )),
  make_option(c("--reseg"), type = "logical",
              default = TRUE,
              help = "Use RRBS ASCAT.m re-segmented (multi-sample) outputs [default %default]"),
  make_option(c("--drop_x"), type = "logical",
              default = FALSE,
              help = "Drop chromosome X from the plot [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!is.null(opt$patient_id) && nzchar(opt$patient_id)) {
  opt$patient_id <- paste(resolve_tx_patient_id(trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE))), .tx_map), collapse = ",")
}

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) {
  stop("Working directory does not exist: ", work_dir)
}
setwd(work_dir)

output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tx_dir <- path.expand(opt$tx_dir)
patient_ids <- trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE)))
sample_ids_opt <- if (is.null(opt$sample_id)) {
  NULL
} else {
  trimws(unlist(strsplit(opt$sample_id, ",", fixed = TRUE)))
}
reseg <- isTRUE(opt$reseg)
drop_x <- isTRUE(opt$drop_x)

################################################################################
# Compare ASCAT.m (RRBS) vs WGS ASCN per patient / region
################################################################################

for (pid_in in patient_ids) {
  ids <- tryCatch(
    resolve_patient_ids(pid_in, tx_dir, opt$clinical_info),
    error = function(e) {
      warning(conditionMessage(e))
      NULL
    }
  )
  if (is.null(ids)) next
  patient_id <- ids$ltx
  cruk_id <- ids$cruk
  path_patient <- file.path(tx_dir, patient_id)
  if (!dir.exists(path_patient)) {
    warning("Patient directory missing, skipping: ", path_patient)
    next
  }

  # Discover regions that have both RRBS and WGS ASCAT outputs unless specified
  if (is.null(sample_ids_opt)) {
    wgs_dirs <- list.dirs(
      file.path(path_patient, "WGS"), full.names = FALSE, recursive = FALSE
    )
    sample_ids <- wgs_dirs[nzchar(wgs_dirs)]
  } else {
    sample_ids <- sample_ids_opt
  }

  for (sample_id in sample_ids) {
    if (reseg) {
      rrbs_path <- file.path(
        path_patient, "Copy_number", sample_id, "re_segmented",
        paste0(patient_id, ".", sample_id, ".ascat.output.resegmented.RData")
      )
    } else {
      rrbs_path <- file.path(
        path_patient, "Copy_number", sample_id, "raw",
        paste0(patient_id, ".", sample_id, ".ascat.output.RData")
      )
      if (!file.exists(rrbs_path)) {
        rrbs_path <- file.path(
          path_patient, "Copy_number", sample_id,
          paste0(patient_id, ".", sample_id, ".ascat.output.RData")
        )
      }
    }
    wgs_path <- file.path(
      path_patient, "WGS", sample_id,
      paste0(sample_id, ".ascat.output.RData")
    )

    if (!file.exists(rrbs_path) || !file.exists(wgs_path)) {
      message(
        "Skipping ", cruk_id, " ", sample_id,
        " (missing RRBS and/or WGS ASCAT output)"
      )
      next
    }

    message("Plotting ASCN comparison: ", cruk_id, " (", patient_id, ") ", sample_id)

    dt_rrbs <- load_ascat_segments(rrbs_path, platform = "RRBS")
    dt_wgs <- load_ascat_segments(wgs_path, platform = "WGS")

    if (drop_x) {
      dt_rrbs <- dt_rrbs[chrom_num != 23]
      dt_wgs <- dt_wgs[chrom_num != 23]
    }

    dt_plot <- build_agreement_plot_dt(dt_rrbs, dt_wgs)
    chr <- chrom_offsets(rbind(dt_rrbs, dt_wgs, fill = TRUE))
    dt_plot <- merge(dt_plot, chr[, .(chrom_num, offset)],
                     by = "chrom_num", all.x = TRUE)

    outfile <- file.path(
      output_dir,
      paste0(cruk_id, "_", sample_id, "_WGS_RRBS_CN_comparison.pdf")
    )
    plot_ascat_vs_wgs_cn(dt_plot, chr, outfile)
    message("Wrote ", outfile)
  }
}

message("Done.")
