#' @description TRACERx RRBS vs WGS allele-specific copy-number intersections
#' for the samples with matched WGS ASCAT outputs (currently n = 7).
#'
#' @author Lili Larose Cadieux
#' @return Writes inputs/TRACERx_RRBS_WGS_overlap.txt (and .RData) for downstream
#' plotting / joint purity-ploidy models.

# Load library
x <- c("data.table", "IRanges", "GenomicRanges", "optparse")
invisible(lapply(x, function(y)
  suppressWarnings(
    suppressMessages(
      library(y,
              character.only = TRUE,
              quietly = TRUE,
              warn.conflicts = FALSE))
    )
  )
)
rm(x)

################################################################################
# Source functions (same logic as rrbs_dnaseq_intersections.R; assay allows WGS)
################################################################################

#' Parse ASCAT output for RRBS or DNA-seq (WES/WGS)
parse_ascat_output <- function(ascat_output_file, patient_id, sample_id, assay = NULL) {
  stopifnot(is.character(assay), is.element(assay, c("wes", "wgs", "rrbs")))
  load(ascat_output_file)

  purity <- ascat.output$aberrantcellfraction[[1]]
  ploidy <- ascat.output$ploidy[[1]]
  fit <- ascat.output$goodnessOfFit[[1]]

  segments_gr <- GenomicRanges::makeGRangesFromDataFrame(
    ascat.output$segments,
    start.field = "startpos",
    end.field = "endpos",
    keep.extra.columns = TRUE
  )

  fraction_loh <- get_loh_fraction(segments_gr)
  wgd <- classify_wgd(ploidy, fraction_loh)

  parsed_data <- list(
    patient_id, sample_id, purity, ploidy, fraction_loh, wgd, fit,
    segments_gr, ascat.output
  )
  names(parsed_data) <- c(
    "patient_id",
    paste0(assay, "_sample"),
    paste0(assay, "_purity"),
    paste0(assay, "_ploidy"),
    paste0(assay, "_loh"),
    paste0(assay, "_wgd"),
    paste0(assay, "_fit"),
    paste0(assay, "_segs"),
    paste0(assay, "_ascat")
  )
  return(parsed_data)
}

get_loh_fraction <- function(segments_gr) {
  segments_gr$loh <- segments_gr$nMajor >= 1 & segments_gr$nMinor == 0
  bases_loh <- sum(width(segments_gr[segments_gr$loh]))
  total_bases <- sum(width(segments_gr))
  return(bases_loh / total_bases)
}

classify_wgd <- function(ploidy, fraction_loh) {
  # PCAWG evolutionary histories WGD classifier
  hom <- fraction_loh
  .classWgd <- function(ploidy, hom) 2.9 - 2 * hom <= ploidy
  .classWgd(ploidy, hom)
}

get_intersection <- function(rrbs_segs, dnaseq_segs) {
  # findOverlapPairs is re-exported via GenomicRanges attach (IRanges)
  intersection <- pintersect(findOverlapPairs(rrbs_segs, dnaseq_segs))
  rrbs_meta <- mergeByOverlaps(rrbs_segs, intersection)
  rrbs_int <- rrbs_meta$intersection
  mcols(rrbs_int) <- mcols(rrbs_meta$rrbs_segs)
  dnaseq_meta <- mergeByOverlaps(dnaseq_segs, intersection)
  dnaseq_int <- dnaseq_meta$intersection
  mcols(dnaseq_int) <- mcols(dnaseq_meta$dnaseq_segs)
  dnaseq_int <- data.frame(dnaseq_int)
  rrbs_int <- data.frame(rrbs_int)
  combined <- merge(rrbs_int, dnaseq_int, by = c("seqnames", "start", "end"))
  stopifnot(all(combined$width.x == combined$width.y))
  return(combined)
}

calculate_copy_overlap <- function(combined, within_1 = FALSE) {
  copy_overlap_bool <-
    (combined$nMajor.x == combined$nMajor.y) &
    (combined$nMinor.x == combined$nMinor.y)
  if (within_1) {
    major_within_1 <-
      (abs(combined$nMajor.x - combined$nMajor.y) <= 1) &
      (combined$nMinor.x == combined$nMinor.y)
    minor_within_1 <-
      (combined$nMajor.x == combined$nMajor.y) &
      (abs(combined$nMinor.x - combined$nMinor.y) <= 1)
    copy_overlap_bool <- major_within_1 | minor_within_1
  }
  bases_copy_overlap <- sum(combined[copy_overlap_bool, ]$width.x)
  bases_total_intersection <- sum(combined$width.x)
  perc_copy_overlap <- bases_copy_overlap * 100 / bases_total_intersection
  list(
    bases_copy_overlap = bases_copy_overlap,
    bases_total_intersection = bases_total_intersection,
    perc_copy_overlap = perc_copy_overlap
  )
}

#' Compare matched RRBS and WGS ASCAT segments for one sample
main_compare_rrbs_wgs <- function(patient_id, rrbs_ascat_file, rrbs_sample_name,
                                  wgs_ascat_file, wgs_sample_name,
                                  cruk_id = NA_character_) {
  rrbs_data <- parse_ascat_output(
    rrbs_ascat_file, patient_id, rrbs_sample_name, assay = "rrbs"
  )
  wgs_data <- parse_ascat_output(
    wgs_ascat_file, patient_id, wgs_sample_name, assay = "wgs"
  )

  # rename WGS segs slot for get_intersection second argument naming clarity
  intersection <- get_intersection(rrbs_data$rrbs_segs, wgs_data$wgs_segs)
  copy_overlap_results <- calculate_copy_overlap(intersection)
  copy_overlap_within_1 <- calculate_copy_overlap(intersection, within_1 = TRUE)

  wgs_segs_bases <- sum(width(wgs_data$wgs_segs))
  rrbs_segs_bases <- sum(width(rrbs_data$rrbs_segs))

  data.table(
    sample_id = paste(ifelse(is.na(cruk_id), patient_id, cruk_id),
                      rrbs_sample_name, sep = "-"),
    patient_id = patient_id,
    tx_id = cruk_id,
    rrbs_sample = rrbs_sample_name,
    wgs_sample = wgs_sample_name,
    purity_a = as.numeric(wgs_data$wgs_purity),
    purity_b = as.numeric(rrbs_data$rrbs_purity),
    ploidy_a = as.numeric(wgs_data$wgs_ploidy),
    ploidy_b = as.numeric(rrbs_data$rrbs_ploidy),
    name_a = "WGS",
    name_b = "RRBS",
    frac_loh_a = as.numeric(wgs_data$wgs_loh),
    frac_loh_b = as.numeric(rrbs_data$rrbs_loh),
    wgs_wgd = as.logical(wgs_data$wgs_wgd),
    rrbs_wgd = as.logical(rrbs_data$rrbs_wgd),
    wgd_concordant = as.logical(rrbs_data$rrbs_wgd == wgs_data$wgs_wgd),
    bases_overlap = as.numeric(copy_overlap_results$bases_copy_overlap),
    bases_intersection = as.numeric(copy_overlap_results$bases_total_intersection),
    perc_overlap = as.numeric(copy_overlap_results$perc_copy_overlap),
    bases_overlap_1 = as.numeric(copy_overlap_within_1$bases_copy_overlap),
    bases_intersection_1 = as.numeric(copy_overlap_within_1$bases_total_intersection),
    perc_overlap_1 = as.numeric(copy_overlap_within_1$perc_copy_overlap),
    wgs_segs_bases = as.numeric(wgs_segs_bases),
    rrbs_segs_bases = as.numeric(rrbs_segs_bases),
    rrbs_ascat_file = rrbs_ascat_file,
    wgs_ascat_file = wgs_ascat_file
  )
}

resolve_rrbs_ascat <- function(txmeth_dir, patient_id, sample_id) {
  reseg <- file.path(
    txmeth_dir, patient_id, "Copy_number", sample_id, "re_segmented",
    paste0(patient_id, ".", sample_id, ".ascat.output.resegmented.RData")
  )
  if (file.exists(reseg)) return(reseg)
  raw <- file.path(
    txmeth_dir, patient_id, "Copy_number", sample_id, "raw",
    paste0(patient_id, ".", sample_id, ".ascat.output.RData")
  )
  if (file.exists(raw)) return(raw)
  NA_character_
}

resolve_wgs_ascat <- function(txmeth_dir, patient_id, sample_id) {
  # same path convention as figure1H-K.R
  main <- file.path(
    txmeth_dir, patient_id, "WGS", sample_id,
    paste0(sample_id, ".ascat.output.RData")
  )
  if (file.exists(main)) return(main)
  NA_character_
}

################################################################################
# parse command line arguments
################################################################################

option_list <- list(
  make_option(c("--wdir"),
              type = "character", default = NULL,
              help = "Working repo base directory (figure1H-K)."),
  make_option(c("--txmeth_dir"),
              type = "character", default = NULL,
              help = "TRACERx methylation (RRBS) + WGS ASCAT repository.")
)
opt <- parse_args(OptionParser(option_list = option_list))

work_dir <- ifelse(is.null(opt$wdir),
                   path.expand("~/data/CAMDAC_paper/figure1H-K"),
                   opt$wdir)
if (!file.exists(work_dir)) {
  stop(paste0("Working directory does not exist: ", work_dir))
}
setwd(work_dir)

txmeth_dir <- ifelse(is.null(opt$txmeth_dir),
                     path.expand("~/rdeocd/txMeth"),
                     opt$txmeth_dir)
if (!file.exists(txmeth_dir)) {
  stop(paste0("Please provide a valid TRACERx methylation path: ", txmeth_dir))
}

in_dir <- "inputs"
out_dir <- "outputs"
if (!dir.exists(in_dir)) dir.create(in_dir)
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# Identify TRACERx samples with WGS ASCAT solutions
################################################################################

f_nm <- file.path(
  txmeth_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
patient_info <- fread(f_nm, header = TRUE, stringsAsFactors = FALSE)
patient_info <- patient_info[!grepl("LN", Sample)]
patient_info <- patient_info[paste(CRUKID, Sample, sep = "_") != "CRUK0084_R4"]

# prefer clinical WGS_rho when present; also keep any sample with a WGS ASCAT file
wgs_candidates <- patient_info[!is.na(WGS_rho),
                               .(Patient, CRUKID, Sample, LTX_ID, WGS_rho, rho_reseg)]

# fall back: scan for WGS folders if clinical annotation is incomplete
if (nrow(wgs_candidates) == 0) {
  warning("No WGS_rho in clinical table; scanning WGS ASCAT files instead.")
}

# ensure every candidate has both ASCAT files; drop those that do not
wgs_candidates[, rrbs_ascat_file :=
  mapply(resolve_rrbs_ascat, txmeth_dir, Patient, Sample)]
wgs_candidates[, wgs_ascat_file :=
  mapply(resolve_wgs_ascat, txmeth_dir, Patient, Sample)]
wgs_candidates <- wgs_candidates[
  !is.na(rrbs_ascat_file) & !is.na(wgs_ascat_file)
]

if (nrow(wgs_candidates) == 0) {
  stop("No TRACERx samples with both RRBS and WGS ASCAT outputs were found.")
}

cat("TRACERx RRBS-WGS pairs to process:\n")
print(wgs_candidates[, .(Patient, CRUKID, Sample, WGS_rho, rho_reseg)])

################################################################################
# Run intersections
################################################################################

rows <- vector("list", nrow(wgs_candidates))
for (i in seq_len(nrow(wgs_candidates))) {
  pid <- wgs_candidates$Patient[i]
  sid <- wgs_candidates$Sample[i]
  cat(sprintf("[%d/%d] %s %s\n", i, nrow(wgs_candidates), pid, sid))
  rows[[i]] <- main_compare_rrbs_wgs(
    patient_id = pid,
    rrbs_ascat_file = wgs_candidates$rrbs_ascat_file[i],
    rrbs_sample_name = sid,
    wgs_ascat_file = wgs_candidates$wgs_ascat_file[i],
    wgs_sample_name = sid,
    cruk_id = wgs_candidates$CRUKID[i]
  )
}

dt_tracerx_wgs <- rbindlist(rows, use.names = TRUE, fill = TRUE)

# CRC-compatible columns first, then TRACERx extras
crc_cols <- c(
  "sample_id", "purity_a", "purity_b", "ploidy_a", "ploidy_b",
  "name_a", "name_b", "frac_loh_a", "frac_loh_b",
  "bases_overlap", "bases_intersection", "perc_overlap",
  "bases_overlap_1", "bases_intersection_1", "perc_overlap_1"
)
extra_cols <- setdiff(names(dt_tracerx_wgs), crc_cols)
dt_out <- dt_tracerx_wgs[, .SD, .SDcols = c(crc_cols, extra_cols)]

out_txt <- file.path(in_dir, "TRACERx_RRBS_WGS_overlap.txt")
out_rdata <- file.path(in_dir, "dt_tracerx_rrbs_wgs.RData")
fwrite(dt_out, out_txt, sep = "\t")
save(dt_tracerx_wgs, file = out_rdata)

cat("\nWrote", nrow(dt_out), "rows to\n ", out_txt, "\n ", out_rdata, "\n")
print(dt_out[, .(sample_id, purity_a, purity_b, ploidy_a, ploidy_b,
                 perc_overlap, perc_overlap_1, wgd_concordant)])

# END
