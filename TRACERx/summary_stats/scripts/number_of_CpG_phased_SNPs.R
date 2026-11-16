#!/usr/bin/env Rscript

#' @description CLI summary of phased SNP/CpG counts and raw allele-count
#' coverage in the CAMDAC dataset:
#' 1. Unique phased SNPs and CpGs.
#' 2. Mean and median CpG coverage overall, outside SNPs, and at SNP-overlapping CpGs.
#' 3. Mean and median total and informative SNP depth.
#'
#' @author Lili Larose Cadieux

# CLI
option_list <- list(
  optparse::make_option(
    "--tx_dir", type = "character", default = "~/rdeocd/txMeth",
    help = "CAMDAC data root [default %default]"
  ),
  optparse::make_option(
    "--work_dir", type = "character",
    default = "~/data/CAMDAC_paper/summary_stats",
    help = "Working directory [default %default]"
  ),
  optparse::make_option(
    "--output_dir", type = "character", default = "outputs",
    help = "Absolute path or path below work_dir [default %default]"
  ),
  optparse::make_option(
    "--patients", type = "character", default = "",
    help = "Optional comma-separated CRUKID (or Patient) IDs"
  ),
  optparse::make_option(
    "--samples", type = "character", default = "",
    help = "Optional comma-separated sample IDs (e.g. N,R1,R2)"
  ),
  optparse::make_option(
    "--schema_only", action = "store_true", default = FALSE,
    help = "Print required allele-count columns and exit"
  ),
  optparse::make_option(
    "--source", type = "character", default = "allelecounts",
    help = paste0(
      "Data source: allelecounts (bulk gr_tumour under Allelecounts/) ",
      "or asm (phased ASM outputs) [default %default]"
    )
  ),
  optparse::make_option(
    "--clinical_info", type = "character", default = "",
    help = "Optional clinical TSV for Patient<->CRUKID map"
  )
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
if (!is.null(opt$patients) && nzchar(opt$patients)) {
  opt$patients <- paste(resolve_tx_patient_id(trimws(unlist(strsplit(opt$patients, ",", fixed = TRUE))), .tx_map), collapse = ",")
}
opt$source <- tolower(trimws(opt$source))
if (!opt$source %in% c("allelecounts", "asm")) {
  stop("--source must be 'allelecounts' or 'asm'")
}

if (opt$schema_only) {
  cat(
    "Bulk CpG coverage: total_counts_m\n",
    "ASM phased CpG coverage: total_counts_m_1 / total_counts_m_2 ",
    "(genotype_1 / genotype_2, not ref/alt labels)\n",
    "SNP depth: total_depth (includes other_counts)\n",
    "Informative SNP depth: total_counts = ref_counts + alt_counts\n",
    "Classification: CpG rows have width > 1; SNP-overlap has non-NA POS\n",
    sep = ""
  )
  quit(save = "no", status = 0)
}

work_dir <- path.expand(opt$work_dir)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
output_dir <- if (grepl("^(/|~)", opt$output_dir)) {
  path.expand(opt$output_dir)
} else {
  file.path(work_dir, opt$output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_tag <- if (opt$source == "asm") "_ASM" else ""
sink(
  file.path(output_dir, paste0("log", output_tag, ".out")),
  append = FALSE, split = TRUE
)

# Load library
x <- c("data.table", "dplyr", "optparse")

invisible(lapply(x, function(y) {
  suppressWarnings(
    suppressMessages(
      library(y, character.only = TRUE, quietly = TRUE,
          warn.conflicts = FALSE)
        )
      )
    })     
  )
rm(x)

# Set base directory
base_dir <- path.expand(opt$tx_dir)

split_arg <- function(value) {
  if (is.null(value) || !nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

requested_patients <- split_arg(opt$patients)
requested_samples <- split_arg(opt$samples)

# Function to get all patient directories
get_patient_dirs <- function(base_dir) {
  patient_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  patient_dirs <- patient_dirs[grepl("^[A-Za-z]+[0-9]", basename(patient_dirs))]
  if (length(requested_patients)) {
    patient_dirs <- patient_dirs[basename(patient_dirs) %in% requested_patients]
  }
  return(patient_dirs)
}

allowed_samples <- c("N", paste0("R", 1:8))

filter_sample_names <- function(sample_names) {
  sample_names <- sample_names[sample_names %in% allowed_samples]
  if (length(requested_samples)) {
    sample_names <- sample_names[sample_names %in% requested_samples]
  }
  sample_names
}

# ASM sample directories under Allele_specific_analysis/
get_asm_sample_dirs <- function(patient_dir) {
  asm_dir <- file.path(patient_dir, "Allele_specific_analysis")
  if (!dir.exists(asm_dir)) return(character())
  sample_names <- filter_sample_names(list.dirs(asm_dir, full.names = FALSE, recursive = FALSE))
  file.path(asm_dir, sample_names)
}

# Bulk methylation allele-count samples under Allelecounts/
get_allelecounts_sample_dirs <- function(patient_dir) {
  ac_dir <- file.path(patient_dir, "Allelecounts")
  if (!dir.exists(ac_dir)) return(character())
  sample_names <- filter_sample_names(list.dirs(ac_dir, full.names = FALSE, recursive = FALSE))
  file.path(ac_dir, sample_names)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

default_clinical_info <- function(tx_dir) {
  file.path(
    tx_dir,
    "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
  )
}

map_patient_to_cruk <- function(dt, tx_dir, clinical_info = "") {
  clinical_info <- if (nzchar(clinical_info)) {
    path.expand(clinical_info)
  } else {
    default_clinical_info(tx_dir)
  }
  if (!file.exists(clinical_info)) {
    return(dt)
  }
  map <- unique(data.table::fread(clinical_info)[, .(ltx_id = Patient, CRUKID)])
  dt[, ltx_id := patient]
  dt <- merge(dt, map, by.x = "patient", by.y = "ltx_id", all.x = TRUE)
  dt[, patient := data.table::fifelse(!is.na(CRUKID), CRUKID, patient)]
  dt[, CRUKID := NULL]
  dt
}

allelecount_coverage_cols <- c(
  "CHR", "start", "end", "width", "POS",
  "total_counts_m", "total_depth", "total_counts"
)

subset_allelecount_columns <- function(data) {
  data <- standardize_allelecount_dt(data)
  missing <- setdiff(allelecount_coverage_cols, names(data))
  if (length(missing)) {
    stop("Missing allele-count columns: ", paste(missing, collapse = ", "))
  }
  data[, ..allelecount_coverage_cols]
}

asm_phased_locus_cols <- c("CHR", "start", "end", "width", "POS")

subset_asm_phased_coverage_columns <- function(data) {
  data <- standardize_allelecount_dt(data)
  missing_locus <- setdiff(asm_phased_locus_cols, names(data))
  if (length(missing_locus)) {
    stop(
      "Missing phased ASM locus columns: ",
      paste(missing_locus, collapse = ", ")
    )
  }
  if (all(c("total_counts_m_1", "total_counts_m_2") %in% names(data))) {
    out <- data[, c(
      ..asm_phased_locus_cols, "total_counts_m_1", "total_counts_m_2"
    ), with = FALSE]
  } else if (all(c("total_counts_m_n_1", "total_counts_m_n_2") %in% names(data))) {
    out <- data[, c(
      ..asm_phased_locus_cols, "total_counts_m_n_1", "total_counts_m_n_2"
    ), with = FALSE]
    data.table::setnames(
      out,
      c("total_counts_m_n_1", "total_counts_m_n_2"),
      c("total_counts_m_1", "total_counts_m_2")
    )
  } else {
    stop(
      "Missing phased ASM coverage columns: need total_counts_m_1/_2 ",
      "or total_counts_m_n_1/_2"
    )
  }
  out
}

empty_phased_genotype_coverage <- function() {
  list(
    mean_cpg_coverage_genotype_1 = NA_real_,
    median_cpg_coverage_genotype_1 = NA_real_,
    mean_cpg_coverage_genotype_2 = NA_real_,
    median_cpg_coverage_genotype_2 = NA_real_,
    mean_cpg_coverage_genotype_1_non_snp = NA_real_,
    median_cpg_coverage_genotype_1_non_snp = NA_real_,
    mean_cpg_coverage_genotype_2_non_snp = NA_real_,
    median_cpg_coverage_genotype_2_non_snp = NA_real_,
    mean_cpg_coverage_genotype_1_snp = NA_real_,
    median_cpg_coverage_genotype_1_snp = NA_real_,
    mean_cpg_coverage_genotype_2_snp = NA_real_,
    median_cpg_coverage_genotype_2_snp = NA_real_
  )
}

analyze_asm_phased_coverage <- function(data) {
  data <- subset_asm_phased_coverage_columns(data)
  data[, phased_cov := data.table::fcoalesce(total_counts_m_1, 0) +
    data.table::fcoalesce(total_counts_m_2, 0)]

  cpgs <- data[
    !is.na(width) & width > 1 & phased_cov > 0,
    .(
      coverage_g1 = safe_mean(total_counts_m_1),
      coverage_g2 = safe_mean(total_counts_m_2),
      snp_overlapping = !all(is.na(POS))
    ),
    by = .(CHR, start, end)
  ]

  list(
    mean_cpg_coverage_genotype_1 = safe_mean(cpgs$coverage_g1),
    median_cpg_coverage_genotype_1 = safe_median(cpgs$coverage_g1),
    mean_cpg_coverage_genotype_2 = safe_mean(cpgs$coverage_g2),
    median_cpg_coverage_genotype_2 = safe_median(cpgs$coverage_g2),
    mean_cpg_coverage_genotype_1_non_snp = safe_mean(
      cpgs[snp_overlapping == FALSE]$coverage_g1
    ),
    median_cpg_coverage_genotype_1_non_snp = safe_median(
      cpgs[snp_overlapping == FALSE]$coverage_g1
    ),
    mean_cpg_coverage_genotype_2_non_snp = safe_mean(
      cpgs[snp_overlapping == FALSE]$coverage_g2
    ),
    median_cpg_coverage_genotype_2_non_snp = safe_median(
      cpgs[snp_overlapping == FALSE]$coverage_g2
    ),
    mean_cpg_coverage_genotype_1_snp = safe_mean(
      cpgs[snp_overlapping == TRUE]$coverage_g1
    ),
    median_cpg_coverage_genotype_1_snp = safe_median(
      cpgs[snp_overlapping == TRUE]$coverage_g1
    ),
    mean_cpg_coverage_genotype_2_snp = safe_mean(
      cpgs[snp_overlapping == TRUE]$coverage_g2
    ),
    median_cpg_coverage_genotype_2_snp = safe_median(
      cpgs[snp_overlapping == TRUE]$coverage_g2
    )
  )
}

# Minimal columns from Allelecounts *.SNPs.CpGs.all.sorted.RData (not Methylation gr_*_m).
sorted_allelecounts_gr_to_dt <- function(gr) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Install GenomicRanges for --source allelecounts")
  }
  mc <- S4Vectors::mcols(gr)
  data.table::data.table(
    CHR = as.character(GenomicRanges::seqnames(gr)),
    start = as.integer(GenomicRanges::start(gr)),
    end = as.integer(GenomicRanges::end(gr)),
    width = as.integer(GenomicRanges::width(gr)),
    POS = mc$POS,
    total_counts_m = mc$total_counts_m,
    total_depth = mc$total_depth,
    total_counts = mc$total_counts
  )
}

load_sorted_allelecounts <- function(path, sample) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  prefer <- if (sample == "N") {
    c("gr_normal", "gr_tumour")
  } else {
    c("gr_tumour", "gr_normal")
  }
  pick <- intersect(prefer, object_names)[1]
  if (is.na(pick) || !nzchar(pick)) {
    gr_hits <- object_names[
      vapply(object_names, function(n) inherits(env[[n]], "GRanges"), logical(1))
    ]
    if (!length(gr_hits)) {
      stop("No GRanges object in Allelecounts file: ", path)
    }
    pick <- gr_hits[1]
  }
  gr <- env[[pick]]
  dt <- sorted_allelecounts_gr_to_dt(gr)
  rm(gr, env)
  invisible(gc())
  list(data = dt, name = pick)
}

load_rdata_table <- function(
    path, preferred = c("gr_tumour", "df_phased", "df", "dt", "allelecounts")) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  object_names <- unique(c(intersect(preferred, object_names), object_names))
  for (object_name in object_names) {
    object <- env[[object_name]]
    if (inherits(object, "GRanges")) {
      return(list(data = object, name = object_name))
    }
    if (is.data.frame(object) || data.table::is.data.table(object)) {
      return(list(data = data.table::as.data.table(object), name = object_name))
    }
  }
  stop("No tabular or GRanges object found in ", path)
}

find_asm_allelecounts_file <- function(sample_dir) {
  patient_dir <- dirname(dirname(sample_dir))
  sample <- basename(sample_dir)
  roots <- c(
    sample_dir,
    file.path(patient_dir, "Methylation", sample),
    file.path(patient_dir, "Copy_number", sample)
  )
  roots <- roots[dir.exists(roots)]
  hits <- unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = "allele.*counts.*[.]RData$",
      full.names = TRUE, recursive = TRUE, ignore.case = TRUE
    )
  }), use.names = FALSE)
  if (!length(hits)) return(NA_character_)
  sort(unique(hits))[1]
}

find_bulk_allelecounts_file <- function(patient_dir, sample) {
  ac_dir <- file.path(patient_dir, "Allelecounts", sample)
  if (!dir.exists(ac_dir)) return(NA_character_)
  patient <- basename(patient_dir)
  preferred <- file.path(
    ac_dir,
    paste0(patient, ".", sample, ".SNPs.CpGs.all.sorted.RData")
  )
  if (file.exists(preferred)) return(preferred)
  hits <- list.files(
    ac_dir,
    pattern = "[.]SNPs[.]CpGs[.]all[.]sorted[.]RData$",
    full.names = TRUE
  )
  if (!length(hits)) return(NA_character_)
  sort(hits)[1]
}

standardize_allelecount_dt <- function(data) {
  data <- data.table::copy(data)
  if (!"CHR" %in% names(data) && "seqnames" %in% names(data)) {
    data[, CHR := as.character(seqnames)]
  }
  if (!"CHR" %in% names(data) && "chrom" %in% names(data)) {
    data[, CHR := as.character(chrom)]
  }
  data
}

analyze_coverage <- function(data) {
  data <- subset_allelecount_columns(data)

  # Collapse first so CpGs linked to multiple SNP records are not overweighted.
  cpgs <- data[
    !is.na(width) & width > 1 & !is.na(total_counts_m),
    .(
      coverage = safe_mean(total_counts_m),
      snp_overlapping = !all(is.na(POS))
    ),
    by = .(CHR, start, end)
  ]
  snps <- data[
    !is.na(POS),
    .(
      depth = safe_mean(total_depth),
      informative_depth = safe_mean(total_counts)
    ),
    by = .(CHR, POS)
  ]

  list(
    mean_cpg_coverage_all = safe_mean(cpgs$coverage),
    median_cpg_coverage_all = safe_median(cpgs$coverage),
    mean_cpg_coverage_non_snp = safe_mean(
      cpgs[snp_overlapping == FALSE]$coverage
    ),
    median_cpg_coverage_non_snp = safe_median(
      cpgs[snp_overlapping == FALSE]$coverage
    ),
    mean_cpg_coverage_snp = safe_mean(
      cpgs[snp_overlapping == TRUE]$coverage
    ),
    median_cpg_coverage_snp = safe_median(
      cpgs[snp_overlapping == TRUE]$coverage
    ),
    mean_snp_depth = safe_mean(snps$depth),
    median_snp_depth = safe_median(snps$depth),
    mean_informative_snp_depth = safe_mean(snps$informative_depth),
    median_informative_snp_depth = safe_median(snps$informative_depth),
    n_cpg_coverage_all = nrow(cpgs),
    n_cpg_coverage_non_snp = cpgs[snp_overlapping == FALSE, .N],
    n_cpg_coverage_snp = cpgs[snp_overlapping == TRUE, .N],
    n_snp_depth = nrow(snps)
  )
}

# Function to load and analyze allele-specific methylation data for a sample
load_and_analyze_asm_data <- function(sample_dir) {
  patient <- basename(dirname(dirname(sample_dir)))
  sample <- basename(sample_dir)
  
  # Look for the RData file
  rdata_file <- file.path(sample_dir,
    paste0("allele_specific_methylation_rates_", sample, ".RData"))
  
  if (!file.exists(rdata_file)) {
    return(NULL)
  }
  
  tryCatch({
    phased <- load_rdata_table(
      rdata_file, preferred = "df_phased_methylation"
    )
    data <- phased$data
    counts <- analyze_unique_counts(data)

    phased_coverage <- empty_phased_genotype_coverage()
    phased_coverage_status <- "ok:df_phased_methylation"
    phased_coverage_result <- tryCatch(
      list(values = analyze_asm_phased_coverage(data)),
      error = function(error) {
        list(error = conditionMessage(error))
      }
    )
    if (!is.null(phased_coverage_result$values)) {
      phased_coverage <- phased_coverage_result$values
    } else {
      phased_coverage_status <- paste0(
        "error:", phased_coverage_result$error
      )
    }

    allelecounts_file <- find_asm_allelecounts_file(sample_dir)
    coverage <- list(
      mean_cpg_coverage_all = NA_real_,
      median_cpg_coverage_all = NA_real_,
      mean_cpg_coverage_non_snp = NA_real_,
      median_cpg_coverage_non_snp = NA_real_,
      mean_cpg_coverage_snp = NA_real_,
      median_cpg_coverage_snp = NA_real_,
      mean_snp_depth = NA_real_,
      median_snp_depth = NA_real_,
      mean_informative_snp_depth = NA_real_,
      median_informative_snp_depth = NA_real_,
      n_cpg_coverage_all = NA_integer_,
      n_cpg_coverage_non_snp = NA_integer_,
      n_cpg_coverage_snp = NA_integer_,
      n_snp_depth = NA_integer_
    )
    coverage_status <- "missing allele-count RData"
    if (!is.na(allelecounts_file)) {
      coverage_result <- tryCatch({
        raw <- load_rdata_table(allelecounts_file)
        object_label <- raw$name
        cov_data <- subset_allelecount_columns(raw$data)
        rm(raw)
        list(
          values = analyze_coverage(cov_data),
          status = paste0("ok:", object_label)
        )
      }, error = function(error) {
        list(
          values = NULL,
          status = paste0("error:", conditionMessage(error))
        )
      })
      if (!is.null(coverage_result$values)) coverage <- coverage_result$values
      coverage_status <- coverage_result$status
    }

    result <- c(list(
      patient = patient,
      sample = sample,
      sample_id = paste0(patient, "_", sample),
      unique_snps = counts$snps,
      unique_cpgs = counts$cpgs,
      total_rows = nrow(data),
      phased_file = rdata_file,
      allelecounts_file = allelecounts_file,
      coverage_status = coverage_status,
      phased_coverage_status = phased_coverage_status
    ), coverage, phased_coverage)
    rm(data)
    invisible(gc())
    return(result)
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

analyze_unique_counts_bulk <- function(data) {
  data <- subset_allelecount_columns(data)
  if (is.null(data) || nrow(data) == 0) {
    return(list(snps = 0, cpgs = 0))
  }
  unique_snps <- uniqueN(data[!is.na(POS), .(CHR, POS)])
  unique_cpgs <- uniqueN(data[
    width > 1 & !is.na(total_counts_m) & total_counts_m > 0,
    .(CHR, start, end)
  ])
  list(snps = unique_snps, cpgs = unique_cpgs)
}

load_and_analyze_bulk_allelecounts <- function(patient_dir, sample_dir) {
  patient <- basename(patient_dir)
  sample <- basename(sample_dir)
  allelecounts_file <- find_bulk_allelecounts_file(patient_dir, sample)
  if (is.na(allelecounts_file)) {
    return(NULL)
  }

  tryCatch({
    raw <- load_sorted_allelecounts(allelecounts_file, sample)
    data <- raw$data
    object_label <- raw$name
    total_rows <- nrow(data)
    counts <- analyze_unique_counts_bulk(data)
    coverage_result <- tryCatch({
      list(
        values = analyze_coverage(data),
        status = paste0("ok:", object_label, ":Allelecounts_sorted")
      )
    }, error = function(error) {
      list(
        values = NULL,
        status = paste0("error:", conditionMessage(error))
      )
    })
    coverage <- list(
      mean_cpg_coverage_all = NA_real_,
      median_cpg_coverage_all = NA_real_,
      mean_cpg_coverage_non_snp = NA_real_,
      median_cpg_coverage_non_snp = NA_real_,
      mean_cpg_coverage_snp = NA_real_,
      median_cpg_coverage_snp = NA_real_,
      mean_snp_depth = NA_real_,
      median_snp_depth = NA_real_,
      mean_informative_snp_depth = NA_real_,
      median_informative_snp_depth = NA_real_,
      n_cpg_coverage_all = NA_integer_,
      n_cpg_coverage_non_snp = NA_integer_,
      n_cpg_coverage_snp = NA_integer_,
      n_snp_depth = NA_integer_
    )
    if (!is.null(coverage_result$values)) {
      coverage <- coverage_result$values
    }

    result <- c(list(
      patient = patient,
      sample = sample,
      sample_id = paste0(patient, "_", sample),
      unique_snps = counts$snps,
      unique_cpgs = counts$cpgs,
      total_rows = total_rows,
      phased_file = NA_character_,
      allelecounts_file = allelecounts_file,
      coverage_status = coverage_result$status
    ), coverage)
    rm(data, raw, counts, coverage_result, coverage)
    invisible(gc())
    result
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

# Function to analyze unique SNPs and CpGs (phased ASM methylation rates)
analyze_unique_counts <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(list(snps = 0, cpgs = 0))
  }
  
  # Count unique SNPs (CHR, POS_genotype)
  unique_snps <- data %>%
    select(CHR, POS_genotype) %>%
    distinct() %>%
    nrow()
  
  # Count unique CpGs (CHR, start, end) where total_counts_m_n_2 > 0 and total_counts_m_n_1 > 0
  unique_cpgs <- data %>%
    filter(total_counts_m_n_2 > 0 & total_counts_m_n_1 > 0 & 
      !is.na(total_counts_m_n_2) & !is.na(total_counts_m_n_1)) %>%
    select(CHR, start, end) %>%
    distinct() %>%
    nrow()
  
  return(list(snps = unique_snps, cpgs = unique_cpgs))
}

# Main execution
cat("Data source:", opt$source, "\n")
if (opt$source == "asm") {
  cat("Loading phased ASM methylation and allele-count coverage...\n")
} else {
  cat("Loading bulk Allelecounts *.SNPs.CpGs.all.sorted.RData ...\n")
}

# Get all patient directories
patient_dirs <- get_patient_dirs(base_dir)
cat("Found", length(patient_dirs), "patient directories\n")

# Initialize results storage as named list
all_results <- list()
total_unique_snps <- 0
total_unique_cpgs <- 0
patient_count <- 0
sample_count <- 0

# Process each patient
for (patient_dir in patient_dirs) {
  patient <- basename(patient_dir)
  cat("\nProcessing patient:", patient, "\n")
  
  sample_dirs <- if (opt$source == "asm") {
    get_asm_sample_dirs(patient_dir)
  } else {
    get_allelecounts_sample_dirs(patient_dir)
  }

  if (!length(sample_dirs)) {
    cat("  No sample directories found\n")
    next
  }

  cat("  Found", length(sample_dirs), "sample directories\n")

  for (sample_dir in sample_dirs) {
    sample <- basename(sample_dir)
    cat("  Loading sample:", sample, "... ")

    result <- if (opt$source == "asm") {
      load_and_analyze_asm_data(sample_dir)
    } else {
      load_and_analyze_bulk_allelecounts(patient_dir, sample_dir)
    }
    
    if (!is.null(result)) {
      # Store result in named list
      all_results[[result$sample_id]] <- result
      
      total_unique_snps <- total_unique_snps + result$unique_snps
      total_unique_cpgs <- total_unique_cpgs + result$unique_cpgs
      sample_count <- sample_count + 1
      cat("SUCCESS (", result$unique_snps, "SNPs,", result$unique_cpgs, "CpGs)\n")
    } else {
      cat("FAILED\n")
    }
  }
  
  patient_count <- patient_count + 1
}

# Summary
cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY:\n")
cat("Patients processed:", patient_count, "\n")
cat("Samples loaded:", sample_count, "\n")
cat("Total unique SNPs across all samples:", total_unique_snps, "\n")
cat("Total unique CpGs across all samples:", total_unique_cpgs, "\n")

# Show detailed results
if (length(all_results) > 0) {
  cat("\nDetailed results by sample:\n")
  cat("Sample_ID\t\tUnique_SNPs\tUnique_CpGs\tTotal_Rows\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  for (sample_id in names(all_results)) {
    result <- all_results[[sample_id]]
    cat(sprintf("%-20s\t%d\t\t%d\t\t%d\n", 
                sample_id, 
                result$unique_snps, 
                result$unique_cpgs, 
                result$total_rows))
  }
}

# Save results
cat("\nSaving results...\n")
# Create summary data frame
results_df <- rbindlist(all_results, fill = TRUE)
results_df[, sample_id := NULL]
setorder(results_df, patient, sample)
results_df <- map_patient_to_cruk(results_df, base_dir, opt$clinical_info)

fwrite(
  results_df,
  file = file.path(output_dir, paste0("unique_snps_cpgs_summary", output_tag, ".csv"))
)
fwrite(
  results_df,
  file = file.path(output_dir, paste0("unique_snps_cpgs_summary", output_tag, ".tsv")),
  sep = "\t"
)

summarize_cohort <- function(dt) {
  data.table::data.table(
    samples = nrow(dt),
    samples_with_coverage = sum(grepl("^ok:", dt$coverage_status)),
    mean_unique_snps = safe_mean(dt$unique_snps),
    mean_unique_cpgs = safe_mean(dt$unique_cpgs),
    mean_cpg_coverage_all = safe_mean(dt$mean_cpg_coverage_all),
    median_cpg_coverage_all = safe_mean(dt$median_cpg_coverage_all),
    mean_cpg_coverage_non_snp = safe_mean(dt$mean_cpg_coverage_non_snp),
    median_cpg_coverage_non_snp = safe_mean(dt$median_cpg_coverage_non_snp),
    mean_cpg_coverage_snp = safe_mean(dt$mean_cpg_coverage_snp),
    median_cpg_coverage_snp = safe_mean(dt$median_cpg_coverage_snp),
    mean_cpg_coverage_genotype_1 = safe_mean(dt$mean_cpg_coverage_genotype_1),
    median_cpg_coverage_genotype_1 = safe_mean(dt$median_cpg_coverage_genotype_1),
    mean_cpg_coverage_genotype_2 = safe_mean(dt$mean_cpg_coverage_genotype_2),
    median_cpg_coverage_genotype_2 = safe_mean(dt$median_cpg_coverage_genotype_2),
    mean_cpg_coverage_genotype_1_snp = safe_mean(dt$mean_cpg_coverage_genotype_1_snp),
    median_cpg_coverage_genotype_1_snp = safe_mean(
      dt$median_cpg_coverage_genotype_1_snp
    ),
    mean_cpg_coverage_genotype_2_snp = safe_mean(dt$mean_cpg_coverage_genotype_2_snp),
    median_cpg_coverage_genotype_2_snp = safe_mean(
      dt$median_cpg_coverage_genotype_2_snp
    ),
    mean_snp_depth = safe_mean(dt$mean_snp_depth),
    median_snp_depth = safe_mean(dt$median_snp_depth),
    mean_informative_snp_depth = safe_mean(dt$mean_informative_snp_depth),
    median_informative_snp_depth = safe_mean(dt$median_informative_snp_depth)
  )
}

cohort_summary <- rbindlist(list(
  cbind(cohort = "all", summarize_cohort(results_df)),
  cbind(cohort = "normal", summarize_cohort(results_df[sample == "N"])),
  cbind(cohort = "tumor", summarize_cohort(results_df[sample != "N"]))
), fill = TRUE)
setcolorder(cohort_summary, c("cohort", setdiff(names(cohort_summary), "cohort")))
fwrite(
  cohort_summary,
  file = file.path(output_dir, paste0("cohort_summary", output_tag, ".tsv")),
  sep = "\t"
)
save(
  all_results, results_df, cohort_summary, opt,
  file = file.path(output_dir, paste0("unique_snps_cpgs_analysis", output_tag, ".RData"))
)

cat("Results saved to:\n")
cat("- unique_snps_cpgs_analysis", output_tag, ".RData (R object)\n", sep = "")
cat("- unique_snps_cpgs_summary", output_tag, ".csv (CSV file)\n", sep = "")
cat("- unique_snps_cpgs_summary", output_tag, ".tsv (TSV file)\n", sep = "")
cat("- cohort_summary", output_tag, ".tsv (cohort means)\n", sep = "")

# Compute the number of SNPs per patient and CpG per patient
results_df %>%
  filter(sample == "N") %>%
  summarise(total_snps = mean(unique_snps),
            total_cpgs = mean(unique_cpgs)) %>%
  print()
print(cohort_summary)

cat("\nAnalysis completed successfully!\n") 

# Close the log file and exit cleanly
sink()
cat("Analysis completed. Check outputs/log.out for detailed results.\n")
q(status = 0)