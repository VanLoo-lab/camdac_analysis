#' @description Investigating the effect of copy number on methylation 
#' rate
#' 
#' @author Lili Larose Cadieux

# Load library
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
x <- c("data.table", "stringr", "ggplot2", "gridExtra", "scales", "optparse")

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

################################################################################
# Source useful functions
################################################################################

#' @description Round numerical values to 'n' digits
#' @param x Numerical vector - values to round
#' @param digits Numerical value - the number of decimal digits to keep
#' @return rounded numerical vector
round2 <- function(x, digits) {
   ifelse(
      as.integer(x * (10 ^ (digits + 1))) %% 10 >= 5, 
      ceiling(x * (10 ^ digits))/(10 ^ digits), 
      floor(x * (10 ^ digits))/(10 ^ digits)
    )
}

#' @description Find the modes of a multi-modal distribution
#' @param x Numerical vector of values from which to estimate the mode
#' @return Vector of modes
find_modes <- function(x) {
  modes <- NULL
  for ( i in 2:(length(x) - 1) ){
    if ( (x[i] > x[i - 1] ) & (x[i] > x[i + 1] ) ) {
      modes <- c(modes, i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}

#' @description Plot the bulk methylation rate separated by copy number
#' @param df_mt Purified CAMDAC tumour methylation rate data.table
#' (inlcludes bulk data and copy number information)
#' @param p Numerical value representing the sample tumour purity
#' @param multip Character value with the allele-specific copy number state of interest
#' written as Major + Minor allele copy number (i.e. "1+0").
#' @param chr Character value with the chromosomes to subset to if required.
#' @param OC Numerical value representing the combined overconversion rate and methylation loss/errors
#' estimated from the matched normal methylation rate density distribution.
#' @param UC Numerical value representing the underconversion rate and methylation errors estimated 
#' from the matched normal methylation rate density distribution.
#' @param reseg Logical value representing whether ASCAT.m was run in single or multi-sample mode.
#' reeg = TRUE means the multi-sample appraoch was used.
#' @param patient_id Patient ID
#' @param sample_id For multi-region sequencing data, set the sample ID. Otherwise, the 
#' default is to set the sample ID as "R1".
#' @param plot Logical value representing wheter the figures should be generated
#' @param out_file_prefix
#' @param hdi_width 
#' @return The modes of a multi-modal distribution
plot_histogram <-
  function(df_mt, p, multip,
           chr = NULL, OC, UC, reseg = TRUE,
           patient_id, sample_id = "R1",
           plot = TRUE, out_file_prefix,
           hdi_width = as.numeric(),
           bin_width = as.numeric()) {
  # set null vars
  b0 <- NA
  t0 <- NA
  b1 <- NA
  t1 <- NA 
  
  # rename multiplicity character
  multi <- multip

  # subset to relevant copy number multiplicity
  df_mt_sub <- df_mt[!is.na(multi) & as.character(multi) == multip,]

  # remove chromosome X - we do not want to complicate the plot with 
  # lots of ASM at this stage
  df_mt_sub <- df_mt_sub[chromosome != "X",]

  # subset to relevant chromosome
  if(!is.null(chr)){
    df_mt_sub <- df_mt_sub[chromosome == chr,]
    f_prefix  <- paste0(f_prefix , "_chr", chr)
  }
  
  # transfor multi charactor to factor
  df_mt_sub$multi <- factor(df_mt_sub$multi, levels = multip)

  # get total and allele-specific copy number
  CN <- mean(df_mt_sub$CG_CN)
  nA <- substr(multip, 1, 1)
  nB <- substr(multip, 3, 3)
  
  # Get normal copy number at heterozygous SNPs
  CN_norm <- mean(ifelse(!is.na(df_mt_sub$BAF), 1, 2))
  # CN norm is 1 out of 2 at CG-destroying/-forming heterozygous SNPs
  
  # get expected peaks in m_b
  # calculate the expected peak for hypermethylated CpGs (mn = 0)
  expected_peak_0 <-
    (OC * p * CN + UC * CN_norm * (1 - p)) / (p * CN + CN_norm * (1 - p))
  expected_peak_0_max <-
    (1 * (p + 0.005) * (CN + 0.5) + 0.2 * CN_norm * (1 - (p + 0.005))) /
    ((p + 0.005) * (CN + 0.5) + CN_norm * (1 - (p + 0.005)))
  expected_peak_0_min <-
    (OC * (p - 0.005) * (CN - 0.5) + 0 * CN_norm * (1 - (p - 0.005))) /
    ((p - 0.005) * (CN - 0.5) + CN_norm * (1 - (p - 0.005)))
   
  # calculate the expected peak for hypomethylated CpGs (mn = 1)
  expected_peak_1 <-
    (UC * p * CN + OC * CN_norm * (1 - p)) /
    (p * CN + CN_norm * (1 - p))
  expected_peak_1_max <-
    (0 * (p + 0.005) * (CN + 0.5) + 0.8 * CN_norm * (1 - (p + 0.005))) /
    ((p + 0.005) * (CN + 0.5) + CN_norm * (1 - (p + 0.005)))
  expected_peak_1_min <- (UC * (p - 0.005) * (CN - 0.5) + 1 * CN_norm * (1 - (p - 0.005))) /
    ((p - 0.005) * (CN - 0.5) + CN_norm * (1 - (p - 0.005)))
  
  # subset the data to loci confidently unmethylated in the normal
  df_mt_0 <- df_mt_sub[df_mt_sub$m_n_high < UC + hdi_width,]
  tmp <- df_mt_0[m_b > expected_peak_0 - 0.05, .SD, .SDcols = "m_b"]
  y0 <-
    3 * as.numeric(
      sort(
        table(
          cut(tmp$m_b, breaks = length(seq(expected_peak_0, 1, 0.02)))
        ),
      decreasing=TRUE)[1]
    )
  
  # get normal methylations density
  d <- density(tmp$m_b[tmp$m_b > 0.01], adjust = 1)
  
  # find local density maxima
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to OC and demythylation with age
  idx <-
    which(mymodes < expected_peak_0_max + 0.2 &
      mymodes > expected_peak_0_min - 0.2)
  b0 <- mymodes[idx]

  # pick the highest density mode is there are mode than 1
  if(length(b0) > 1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    b0 <- mymodes[which.max(d$y[mymodes_indices])]  
  }

  # set to NA if there are no suitable modes.
  if(length(b0) < 1) b0 <- NA
  rm(tmp, idx)
  
  # subset the data to loci confidently methylated in the normal
  df_mt_1 <- df_mt_sub[df_mt_sub$m_n_low > OC - hdi_width,]
  tmp <- df_mt_1[m_b < expected_peak_1 + 0.05, .SD, .SDcols = "m_b"]
  y1 <-
    3 * as.numeric(
      sort(
        table(
          cut(tmp$m_b, breaks = length(seq(0, expected_peak_1, 0.02)))
        ),
      decreasing=TRUE)[1]
    )
  
  # get normal methylations density
  d <- density(tmp$m_b[tmp$m_b < 0.98], adjust = 0.8)
  
  # find modes
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to UC
  idx <-
    which(mymodes > expected_peak_1_max - 0.2 &
      mymodes < expected_peak_1_min + 0.2)
  b1 <- mymodes[idx]

  # pick the highest density mode is there are mode than 1
  if(length(b1) > 1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    b1 <- mymodes[which.max(d$y[mymodes_indices])]  
  }

  # set to NA if there are no suitable modes.
  if(length(b1) < 1) b1 <- NA
  rm(tmp, d, mymodes, mymodes_indices)

  # set the colors for plotting
  colss <-
    ifelse(multip %in% c("2+0","4+0"), "dodgerblue",
    ifelse(multip %in% c("2+1","4+2"), "orange3",
    ifelse(multip %in% c("1+1","2+2"), "firebrick", 
    ifelse(multip == "1+0", "grey65", "forestgreen"))))

  if (plot == TRUE) {
    p0 <- ggplot(df_mt_0) +
      geom_histogram(aes(x=m_b, y=..count..),
        fill=colss, color=NA, binwidth=bin_width, alpha=0.5) +
      theme_classic(base_size = 7) +
      coord_cartesian(ylim=c(0,y0),xlim=c(0,1))+
      geom_vline(aes(xintercept = expected_peak_0), size = 0.4) +
      geom_vline(aes(xintercept = expected_peak_0_max),
        linetype = "dashed", size = 0.25) +
      geom_vline(aes(xintercept = expected_peak_0_min),
        linetype = "dashed", size = 0.25)
    ggsave(file = paste0(out_file_prefix, "_bulk_methylation_hyper_multi_",
            nA, "_", nB, ".pdf"),
      plot = p0, device = "pdf", width = 1.5, height = 1, units = "in")
  
    p1 <- ggplot(df_mt_1) +
      geom_histogram(aes(x=m_b, y=..count..),
        fill = colss, color = NA, binwidth = bin_width, alpha = 0.5) +
      theme_classic(base_size = 7) +
      coord_cartesian(ylim = c(0,y1),xlim = c(0,1))+
      geom_vline(aes(xintercept=expected_peak_1), size = 0.4) +
      geom_vline(aes(xintercept=expected_peak_1_max),
        linetype="dashed", size = 0.25) +
      geom_vline(aes(xintercept=expected_peak_1_min),
        linetype="dashed", size = 0.25)
    ggsave(file = paste0(out_file_prefix, "_bulk_methylation_hypo_multi_",
      nA, "_", nB,".pdf"),
      plot = p1, device = "pdf", width = 1.5, height = 1, units = "in")
  }
    
  return(c(b0=b0, expected_b0=expected_peak_0, b1=b1, expected_b1=expected_peak_1))
}

#' @description Estimate over-/under-conversion from a normal methylation vector
#' @param m_n Numerical vector of normal methylation rates
#' @param uc_min Lower bound used to select the under-conversion mode
#' @param uc_max Upper bound used to select the under-conversion mode
#' @param oc_min Lower bound used to select the over-conversion mode
#' @param oc_max Upper bound used to select the over-conversion mode
#' @return Named list with OC and UC
estimate_conversion_rates <- function(m_n,
                                      uc_min = 0, uc_max = 0.04,
                                      oc_min = 0.85, oc_max = 0.98) {
  d <- density(m_n, adjust = 0.95)
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  list(
    OC = mymodes[which(mymodes > oc_min & mymodes < oc_max)],
    UC = mymodes[which(mymodes > uc_min & mymodes < uc_max)]
  )
}

#' @description Plot bulk methylation histograms for selected multiplicities
#' @param df_mt Methylation data.table
#' @param p Sample purity
#' @param top_multi Character vector of multiplicities to plot
#' @param OC Over-conversion rate
#' @param UC Under-conversion rate
#' @param reseg Logical ASCAT.m re-segmentation flag passed to plot_histogram
#' @param patient_id Patient ID
#' @param sample_id Sample ID
#' @param out_dir Output directory for PDFs
#' @param hdi_width HDI width passed to plot_histogram
#' @param bin_width Histogram bin width
#' @param verbose If TRUE, print multiplicity and plot_histogram return value
plot_multiplicities <- function(df_mt, p, top_multi, OC, UC, reseg,
                                patient_id, sample_id, out_dir,
                                hdi_width, bin_width, verbose = FALSE) {
  f_prefix <- file.path(out_dir, paste0(patient_id, "_", sample_id))
  for (multi in top_multi[!is.na(top_multi)]) {
    if (verbose) print(multi)
    cat("Generateing histograms - ")
    x <- plot_histogram(
      df_mt = df_mt, p = p,
      multip = multi, chr = NULL,
      OC = OC, UC = UC, reseg = reseg,
      patient_id = patient_id, sample_id = sample_id,
      plot = TRUE, out_file_prefix = f_prefix,
      hdi_width = hdi_width,
      bin_width = bin_width
    )
    if (verbose) print(x)
    cat("histograms generated\n")
  }
}

#' @description Resolve a path that may be absolute, ~-relative, or relative
resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) {
    return(path)
  }
  file.path(base_dir, path)
}


################################################################################
# Parse command line arguments
################################################################################

option_list <- list(
  make_option(c("--project"), type = "character",
              default = "tracerx",
              help = "Project to run: tracerx or mcrpc [default %default]"),
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figure1B-D_S1G-H",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory (relative to --work_dir or absolute) [default %default]"),
  make_option(c("--patient_id"), type = "character",
              default = NULL,
              help = paste0(
                "Comma-separated patient ID(s). ",
                "Defaults: CRUK0062 (tracerx); DTB-097,DTB-170 (mcrpc)"
              )),
  make_option(c("--tx_dir"), type = "character",
              default = "~/rdeocd/txMeth",
              help = "TRACERx methylation parent directory [default %default]"),
  make_option(c("--reseg"), type = "logical",
              default = TRUE,
              help = "TRACERx: ASCAT.m multi-sample (re-segmented) mode [default %default]"),
  make_option(c("--clinical_info"), type = "character",
              default = paste0(
                "tx100clinicalinfo_RRBS_published_status_",
                "immune_clusters_purity_ploidy_CpG_cov.txt"
              ),
              help = paste0(
                "TRACERx clinical info file within --tx_dir, or absolute path ",
                "[default %default]"
              )),
  make_option(c("--input_dir"), type = "character",
              default = "inputs",
              help = "mCRPC input directory (relative to --work_dir or absolute) [default %default]"),
  make_option(c("--panel"), type = "character",
              default = "prostate_panel.csv.gz",
              help = "mCRPC normal panel file within --input_dir [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!is.null(opt$patient_id) && nzchar(opt$patient_id)) {
  opt$patient_id <- paste(resolve_tx_patient_id(trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE))), .tx_map), collapse = ",")
}

project <- tolower(opt$project)
if (!project %in% c("tracerx", "mcrpc")) {
  stop("--project must be 'tracerx' or 'mcrpc'")
}

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) {
  stop(paste0("Working directory does not exist: ", work_dir))
}
setwd(work_dir)

output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

patient_ids <- if (is.null(opt$patient_id)) {
  if (project == "tracerx") "CRUK0062" else c("DTB-097", "DTB-170")
} else {
  trimws(unlist(strsplit(opt$patient_id, ",", fixed = TRUE)))
}


################################################################################
# Get bulk methylation plots with CN and purity and DMPs
################################################################################

if (project == "tracerx") {
  tx_dir <- path.expand(opt$tx_dir)
  reseg <- opt$reseg
  clinical_info <- resolve_path(opt$clinical_info, tx_dir)

  patient_info <- fread(clinical_info, header = TRUE, stringsAsFactors = FALSE)
  # remove lymph nodes
  patient_info <- patient_info[!grepl("LN", Sample),]

  for (patient_id in patient_ids) {
    path_patient <- file.path(tx_dir, patient_id)
    sample_ids <- patient_info[Patient == patient_id, Sample]

    load(file.path(path_patient, "Methylation/N/gr_normal_m.RData"))
    conv <- estimate_conversion_rates(gr_normal_m$m_n, uc_min = 0)
    rm(gr_normal_m)
    OC <- conv$OC
    UC <- conv$UC

    for (sample_id in sample_ids) {
      # Load the CN and normal annotated tumour data
      if (reseg == FALSE | patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))) {
        load(
          paste0(path_patient, "/Methylation/", sample_id,
                 "/purified_tumour_DMP_calls.RData")
        )
      }
      if (reseg == TRUE & !patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))) {
        load(
          paste0(path_patient, "/Methylation/", sample_id,
                 "/re_segmented/purified_tumour_DMP_calls.RData")
        )
      }

      df_mt <- data.table(df_mt)

      # get multiplicitity of copy number segments outside X
      df_top_multi <-
        data.table(
          sort(df_mt[chromosome != "X", table(multi)],
            decreasing = TRUE))
      colnames(df_top_multi) <- c("multip", "counts")

      # get DMPs at each multiplicitities
      tmp <-
        df_mt[, .(hyper_counts = sum(!is.na(DMP_t) &
                  DMP_t == "hyper" & m_n_high < 0.2),
                 hypo_counts = sum(!is.na(DMP_t) &
                  DMP_t == "hypo" & m_n_low > 0.8)),
              by = .(multi)]

      # annotate DMP counts
      df_top_multi[, c("hyper_counts", "hypo_counts") :=
        tmp[match(df_top_multi$multip, tmp$multi), 2:3]
      ]
      # keep common states (we will be looking at CpGs which are fully 0 or 1 in m_n
      # and looking at DM so need high initial numbers)
      df_top_multi <-
        df_top_multi[counts > 1E5 & hypo_counts > 500 & hypo_counts > 500, ]

      #  add results columns
      df_top_multi[, c("b0", "expected_b0", "b1", "expected_b1") := 0]

      top_multi <- sort(df_top_multi$multip)
      cat("Top CpG copy number states :", sep = "\n")
      cat(top_multi, sep = "\n")

      # Load the copy number segments file
      if (reseg == FALSE | patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))) {
        copy_number_filename <-
          paste(patient_id, sample_id, "ascat.output.RData", sep = ".")
        file_CN <-
          paste0(path_patient, "/Copy_number/", sample_id, "/",
            copy_number_filename)
      }
      if (reseg == TRUE & !patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))) {
        copy_number_filename <-
          paste(patient_id, sample_id, "ascat.output.resegmented.RData",
            sep = ".")
        file_CN <-
          paste0(path_patient, "/Copy_number/", sample_id,
            "/re_segmented/", copy_number_filename)
      }
      load(file_CN)
      p <- as.numeric(ascat.output$aberrantcellfraction)
      rm(ascat.output)

      plot_multiplicities(
        df_mt = df_mt, p = p, top_multi = top_multi,
        OC = OC, UC = UC, reseg = reseg,
        patient_id = patient_id, sample_id = sample_id,
        out_dir = output_dir,
        hdi_width = 0.2, bin_width = 0.02,
        verbose = FALSE
      )
    }
  }

} else if (project == "mcrpc") {
  input_dir <- resolve_path(opt$input_dir, work_dir)
  panel_file <- resolve_path(opt$panel, input_dir)

  for (patient_id in patient_ids) {
    sample_id <- "R1"

    # read in the tumour methylation file
    f_nm <- file.path(input_dir, paste0(patient_id, ".T.pure.csv.gz"))
    wgbs_data <- fread(f_nm)

    # read in methylation panel
    wgbs_normal_panel <- fread(panel_file)

    # rename and subset columns
    colnames(wgbs_normal_panel) <-
      gsub("_n", "_n_o", colnames(wgbs_normal_panel))

    # annotate m_n_low and m_n_high
    wgbs_data <-
      merge(wgbs_data, wgbs_normal_panel,
            by = c("chrom", "start", "end"))
    wgbs_data[, "m_n_low" := m_n_o_low]
    wgbs_data[, "m_n_high" := m_n_o_high]

    # reformat column names
    colnames(wgbs_data)[1] <- "chromosome"
    wgbs_data[, "m_b" := m]

    conv <- estimate_conversion_rates(wgbs_data$m_n, uc_min = 0.001)
    OC <- conv$OC
    UC <- conv$UC

    # get top copy number states
    wgbs_data[, "multi" := paste0(nA, "+", nB)]
    top_multi <-
      names(sort(wgbs_data[, table(multi)], decreasing = TRUE)[1:5])
    cat("Top CpG copy number states :", sep = "\n")
    cat(top_multi, sep = "\n")

    # Remove CpGs that are heterozygous across the infiltrating normals
    wgbs_data <- wgbs_data[CG_CN_n == 2]

    # Extract purity
    p <- as.numeric(unique(wgbs_data$p))

    plot_multiplicities(
      df_mt = wgbs_data, p = p, top_multi = top_multi,
      OC = OC, UC = UC, reseg = FALSE,
      patient_id = patient_id, sample_id = sample_id,
      out_dir = output_dir,
      hdi_width = 1 / 3, bin_width = 0.03,
      verbose = TRUE
    )
  }
}

# END