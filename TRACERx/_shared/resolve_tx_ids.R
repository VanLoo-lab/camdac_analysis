#' TRACERx patient ID helpers: prefer CRUKID in scripts/docs; resolve to
#' internal Patient (LTX*) via the tx100 clinical table when file paths need it.
#'
#' @param clinical_file Path to tx100 clinical table with Patient + CRUKID columns.
#'   Default: ~/rdeocd/txMeth/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt

TX_CLINICAL_DEFAULT <- path.expand(file.path(
  "~", "rdeocd", "txMeth",
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
))

#' Single-region TRACERx patients (no multi-sample resegmentation)
TX_SINGLE_REGION_CRUK <- c("CRUK0010", "CRUK0047")

#' Load Patient <-> CRUKID map from the clinical table
#' @return data.table with unique Patient, CRUKID
load_tx_patient_map <- function(clinical_file = TX_CLINICAL_DEFAULT) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required")
  }
  if (!file.exists(clinical_file)) {
    stop("Clinical table not found: ", clinical_file)
  }
  dt <- data.table::fread(clinical_file, header = TRUE, stringsAsFactors = FALSE)
  if (!all(c("Patient", "CRUKID") %in% names(dt))) {
    stop("Clinical table must contain Patient and CRUKID columns")
  }
  unique(dt[, .(Patient, CRUKID)])
}

#' Resolve one or more CRUKIDs (or already-Patient IDs) to Patient (LTX*) IDs
#' @param ids Character vector of CRUKIDs and/or Patient IDs
#' @param patient_map Optional map from load_tx_patient_map(); loaded if NULL
#' @return Character vector of Patient IDs (same length as ids when 1:1)
resolve_tx_patient_id <- function(ids, patient_map = NULL) {
  ids <- as.character(ids)
  if (is.null(patient_map)) patient_map <- load_tx_patient_map()
  out <- character(length(ids))
  for (i in seq_along(ids)) {
    id <- ids[i]
    if (id %in% patient_map$Patient) {
      out[i] <- id
      next
    }
    hits <- unique(patient_map[CRUKID == id, Patient])
    if (length(hits) == 0) {
      stop("Unknown TRACERx id (not Patient or CRUKID): ", id)
    }
    if (length(hits) > 1) {
      stop("CRUKID ", id, " maps to multiple Patient IDs: ",
           paste(hits, collapse = ", "))
    }
    out[i] <- hits
  }
  out
}

#' Patient IDs for the single-region (no-reseg) TRACERx cases
#' @return Character vector of Patient IDs
tx_single_region_patients <- function(patient_map = NULL) {
  resolve_tx_patient_id(TX_SINGLE_REGION_CRUK, patient_map = patient_map)
}

#' Whether a Patient ID is a single-region (no-reseg) case
is_tx_single_region <- function(patient_id, patient_map = NULL) {
  patient_id %in% tx_single_region_patients(patient_map = patient_map)
}
