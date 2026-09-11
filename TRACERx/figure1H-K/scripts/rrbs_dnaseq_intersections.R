#' @description Comparison of genomic and methylation-derived allele-specific
#' copy numbers across TRACERx
#' 
#' @author Nana Mensah, Lili Larose Cadieux
#' @return Writes all intersection dataframes as RData files to the current directory.
#' These are used downstream as plotting inputs

# Load library
x <- c("data.table", "ggplot2", "gridExtra", "GenomicRanges", "optparse")
invisible(lapply(x, function(y) 
  suppressWarnings(
    suppressMessages(
      library(y,
      character.only = TRUE,
      quietly=TRUE,
      warn.conflicts = FALSE))
      )
    )
  )
rm(x)

################################################################################
# Source functions
################################################################################

#' @title main_compare_ascat_segments
#' Main function to process ascat output files from matched rrbs and wes samples
#' @return a dataframe with one row containing the following fields:
#'  patient_id, rrbs_sample, rrbs_purity, rrbs_ploidy, rrbs_loh, rrbs_wgd, rrbs_fit,
#'  wes_sample, wes_purity, wes_ploidy, wes_loh, wes_wgd, wes_fit, perc_copy_overlap, 
#'  wgd_concordant
#'
#' @param patient_id patient identifier
#' @param rrbs_ascat_file Path to ascat output file for RRBS sample (CAMDAC)
#' @param rrbs_sample_name Name of rrbs sample
#' @param wes_ascat_file Path to ascat output file for WES sample
#' @param wes_sample_name Name of wes sample
main_compare_ascat_segments <-
  function(patient_id, rrbs_ascat_file, rrbs_sample_name, wes_ascat_file, wes_sample_name) {

  rrbs_data <- parse_ascat_output(rrbs_ascat_file, patient_id, rrbs_sample_name, assay="rrbs")
  wes_data <- parse_ascat_output(wes_ascat_file, patient_id, wes_sample_name, assay="wes")
  intersection <- get_intersection(rrbs_data$rrbs_segs, wes_data$wes_segs)
  copy_overlap_results <- calculate_copy_overlap(intersection)
  copy_overlap_within_1 <- calculate_copy_overlap(intersection, within_1=TRUE)
  wgd_concordant <- rrbs_data$rrbs_wgd == wes_data$wes_wgd
  wes_segs_bases=sum(width(wes_data$wes_segs))
  rrbs_segs_bases=sum(width(rrbs_data$rrbs_segs))
  

  results_object <- c(
    within(rrbs_data, rm(rrbs_segs, rrbs_ascat)), # remove large segments and ascat objects
    within(wes_data, rm(wes_segs, wes_ascat, patient_id)), # as above, with patient_id removed to avoid duplicate
    wgd_concordant=wgd_concordant,
    perc_copy_overlap=copy_overlap_results$perc_copy_overlap,
    bases_copy_overlap=copy_overlap_results$bases_copy_overlap,
    bases_intersection=copy_overlap_results$bases_total_intersection,
    wes_segs_bases=wes_segs_bases,
    wes_perc_intersect=copy_overlap_results$bases_total_intersection*100/wes_segs_bases,
    rrbs_segs_bases=rrbs_segs_bases,
    rrbs_perc_intersect=copy_overlap_results$bases_total_intersection*100/rrbs_segs_bases,
    perc_copy_within_1=copy_overlap_within_1$perc_copy_overlap
  )
  
  return(data.frame(results_object))
}

#' Read ascat.output file for a given RRBS or WES sample
#' and extract purity, ploidy, fraction_loh, purity, ploidy and ascat segments  
parse_ascat_output <- function(ascat_output_file, patient_id, sample_id, assay = NULL){
  stopifnot(is.character(assay), is.element(assay, c("wes", "rrbs")))
  load(ascat_output_file) # Brings object ascat.output into environment
  
  # Get variables from ascat output
  purity <- ascat.output$aberrantcellfraction[[1]]
  ploidy <- ascat.output$ploidy[[1]]
  fit <- ascat.output$goodnessOfFit[[1]]
  
  # Convert segments to genomic ranges
  segments_gr <- GenomicRanges::makeGRangesFromDataFrame(
    ascat.output$segments, start.field="startpos", end.field="endpos", keep.extra.columns=TRUE
  )
  
  # Calculate loh fraction and wgd status
  fraction_loh <- get_loh_fraction(segments_gr)
  wgd <- classify_wgd(ploidy, fraction_loh)
  
  # Return parsed data object
  parsed_data <- list(
    patient_id, sample_id, purity, ploidy, fraction_loh, wgd, fit, segments_gr, ascat.output
  )
  names(parsed_data) <- c(
    "patient_id", paste0(assay,"_sample"), paste0(assay,"_purity"), paste0(assay,"_ploidy"), paste0(assay,"_loh"),
    paste0(assay,"_wgd"), paste0(assay,"_fit"), paste0(assay,"_segs"), paste0(assay,"_ascat")
  )
  return(parsed_data)
}

get_loh_fraction <- function(segments_gr){
  segments_gr$loh = segments_gr$nMajor >=1 & segments_gr$nMinor == 0
  bases_loh <- sum(width(segments_gr[segments_gr$loh]))
  total_bases <- sum(width(segments_gr))
  return(bases_loh / total_bases)
}

classify_wgd <- function(ploidy, fraction_loh){
  # Returns TRUE or FALSE
  # WGD status is determined by ploidy and LOH status. We use the classification from PCAWG
  #   evolutionary histories paper, defined here:
  #  https://github.com/gerstung-lab/PCAWG-11/blob/6fc785dafbb54b00111899b40d93221785aff0a7/code/PCAWG-functions.R#L612
  #  In their code, `hom` in PCAWG comes from averageHom function, which is fraction_loh (minor copy-no = 0)
  hom <- fraction_loh
  .classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy
  wgd_status <- .classWgd(ploidy, hom)
  return( wgd_status )
}

get_intersection <- function(rrbs_segs, wes_segs){
  # Get the overlapping copy number segments (intersect) and annotate with data from both.
  #   Segments must be granges objects as formatted by parse_ascat_output()
  #   Returns a list of bases with copy overlapping, bases in intersection and percentage copy overlap

  ## Pairs+pintersect required to break up ranges with multiple overlaps
  intersection <- GenomicRanges::pintersect(findOverlapPairs(rrbs_segs, wes_segs))
  rrbs_meta <- mergeByOverlaps(rrbs_segs, intersection)
  rrbs_int <- rrbs_meta$intersection
  mcols(rrbs_int) <- mcols(rrbs_meta$rrbs_segs)
  wes_meta <- mergeByOverlaps(wes_segs, intersection)
  wes_int <- wes_meta$intersection
  mcols(wes_int) <- mcols(wes_meta$wes_segs)
  wes_int <- data.frame(wes_int)
  rrbs_int <- data.frame(rrbs_int)
  combined <- merge(rrbs_int, wes_int, by=c("seqnames","start","end"))
  stopifnot(all(combined$width.x == combined$width.y)) # Assert base widths from overlap are the same
  return(combined)
}

calculate_copy_overlap <- function(combined, within_1=FALSE){
  # Return the percentage of bases in rrbs-wes segment intersection with
  #  the same copy number state across both assays.
  # nMajor.x = rrbs, nMajor.y = wes. Same for nMinors.
  copy_overlap_bool <- (combined$nMajor.x == combined$nMajor.y) & (combined$nMinor.x == combined$nMinor.y)
  major_bool <- (combined$nMajor.x == combined$nMajor.y)
  # If within_1 set, limit to regions where copy overlap is +/- 1 on either major or minor only
  if (within_1){
    major_within_1 <- (abs(combined$nMajor.x - combined$nMajor.y) <= 1) & (combined$nMinor.x == combined$nMinor.y)
    minor_within_1 <- (combined$nMajor.x == combined$nMajor.y) & (abs(combined$nMinor.x - combined$nMinor.y) <= 1)
    copy_overlap_bool <- major_within_1 | minor_within_1
  }

  # Limit to regions of overlap
  bases_copy_overlap <- sum(combined[copy_overlap_bool,]$width.x)
  bases_total_intersection <- sum(combined$width.x)
  perc_copy_overlap <- bases_copy_overlap * 100 / bases_total_intersection
  return(list(bases_copy_overlap=bases_copy_overlap, bases_total_intersection=bases_total_intersection, perc_copy_overlap=perc_copy_overlap))
}

get_and_write_intersection <- function(patient_id, rrbs_ascat_file, rrbs_sample_name, wes_ascat_file, wes_sample_name){

  # Get intersection
  rrbs_data <- parse_ascat_output(rrbs_ascat_file, patient_id, rrbs_sample_name, assay = "rrbs")
  wes_data <- parse_ascat_output(wes_ascat_file, patient_id, wes_sample_name, assay = "wes")
  intersection <- get_intersection(rrbs_data$rrbs_segs, wes_data$wes_segs)

  # Write output
  save(intersection, file = paste0(rrbs_sample_name, ".RData"))
}

################################################################################
# Run main
################################################################################

# parse options list
option_list <- list(
  make_option(c("--wdir"),
              type = "character", default = NULL,
              help = "Working repo base directory."),
  make_option(c("--txmeth_dir"),
              type = "character", default = NULL,
              help = "TRACERx methylation (RRBS) data repository."),
  make_option(c("--txwes_dir"),
              type = "character", default = NULL,
              help = "TRACERx WES data repository.")
)
opt <- parse_args(OptionParser(option_list = option_list))

# get working directory
work_dir <- ifelse(is.null(opt$wdir), getwd(), opt$wdir)
if (!file.exists(work_dir)) {
  stop(paste0("Working directory does not exist:", work_dir))
}
setwd(work_dir)

# get TRACERx methylation data directory
txmeth_dir <- opt$txmeth_dir
if (is.null(opt$txmeth_dir) | !file.exists(txmeth_dir)) {
  stop(paste0("Please provide a valid path for TRACERx methylation data ", txmeth_dir))
}

# get TRACERx methylation data directory
txwes_dir <- opt$txwes_dir
if (is.null(opt$txwes_dir) | !file.exists(txwes_dir)) {
  stop(paste0("Please provide a valid path for TRACERx WES data ", txwes_dir))
}

 # set input and output directories
in_dir <- "inputs"
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir)

if (interactive()) {
  # set working directory
  work_dir <- "~/data/CAMDAC_paper/figure2"
  setwd(work_dir)

  # set parent directory
  txmeth_dir <- "~/rdeocd/txMeth"
  txwes_dir <- "/path/to/wes"
}

# set vars
reseg <- TRUE
f_nm <-
  file.path(txmeth_dir, 
    "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
  )

# load patient_info
patient_info <- fread(f_nm, header = TRUE, stringsAsFactors = FALSE)
rm(f_nm)

# get patient ids
patient_ids <- unique(patient_info$Patient)

# iterate over patient IDs
for(patient_id in patient_ids) {
  # set the patient path
  path_patient_txmeth <- file.path(txmeth_dir, patient_id, "/")
  path_patient_txwes <- file.path(txwes_dir, patient_id, "/")

  # iterate over sample IDs
  sample_ids <- patient_info[Patient == patient_id, Sample]

  for(sample_id in sample_ids) {
    # Get resegmented sample RRBS-derived ASCAT.m output file
    file_nm <- paste0(patient_id, ".", sample_id, ".ascat.output.resegmented.RData")
    file_pth_txmeth <- 
        file.path(path_patient_txmeth, "Copy_number", sample_id, "re_segmented", file_nm)
      
    if(!file.exists(file_pth_txmeth)){ 
      # Get single sample RRBS-derived ASCAT.m output file
      file_nm <- paste0(patient_id, ".", sample_id, ".ascat.output.RData")
      file_pth_txmeth <- 
        file.path(path_patient_txmeth, "Copy_number", sample_id, "raw", file_nm)
    }
    
    # Get WES-derived ASCAT output file 
    file_nm <- paste0(patient_id, ".", sample_id, ".ascat.output.RData")
    file_pth <-
        file.path(path_patient_txmeth, "Copy_number", sample_id, "re_segmented", file_nm)


    # set input data
    input_data <-
      c("patient_id" = patient_id,
        "rrbs_ascat_file" = file_pth_txmeth,
        "rrbs_sample_name" = sample_id,
        "wes_ascat_file" = "to_do",
        "wes_sample_name" = "to_do")
    
    # run script
    do.call(get_and_write_intersection, input_data[i, ])
  }
}

# END