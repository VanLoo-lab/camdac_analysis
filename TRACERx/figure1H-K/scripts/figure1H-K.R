#! usr/bin/Rscript
 
#' @description Comparison of genomic and methylation-derived allele-specific
#' copy numbers and purity values across TRACERx
#' 
#' @author Lili Larose Cadieux

# Load library
x <- c("data.table", "ggplot2", "gridExtra", "GenomicRanges", "dplyr",
  "optparse")
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

###############################################################################
# parse command line arguments
###############################################################################

if (interactive() == TRUE) {
  # parse options list
  option_list <- list(
    make_option(c("--wdir"),
                type = "character", default = NULL,
                help = "Working repo base directory.")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  # read in the work directory
  work_dir <- ifelse(is.null(opt$wdir), getwd(), opt$wdir)
} else {
  # set working directory
  work_dir <- "~/data/CAMDAC_paper/figure1H-K"
}

# get working directory
if (!file.exists(work_dir)) {
  stop(paste0("Working directory does not exist:", work_dir))
}
setwd(work_dir)

# set input and output directories
in_dir <- "inputs"
out_dir <- "outputs"

# set script directory
tracerx_dir <- "~/rdeocd/txMeth"

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

################################################################################
# Annotate the raw and re-segmented ASCAT.m purity and ploidy values, the
# fraction of the genome with LOH and whole genome doubling status
################################################################################

# load patient_info
f_nm <- paste0("tx100clinicalinfo_RRBS_published_status",
               "_immune_clusters_purity_ploidy_CpG_cov.txt")
patient_info_tmp <- fread(file.path(tracerx_dir, f_nm))
patient_info <- patient_info_tmp[, .SD, .SDcols = 1:11]
rm(f_nm)

# create results vector
patient_info[, c("rho_raw", "psi_grid_raw", "psi_plot_raw", "frac_LOH_raw",
                 "rho_reseg", "psi_grid_reseg", "psi_plot_reseg", "frac_LOH_reseg",
                 "WGS_rho", "WGS_psi_grid", "WGS_psi_plot", "WGS_frac_LOH") :=
               as.numeric(NA)]
patient_info[, c("WGD_raw", "WGD_reseg", "WGS_WGD") := as.logical(NA)]

# get patient ids
patient_files <- list.files(tracerx_dir, recursive = FALSE)
patient_ids <- grep("LTX", patient_files, value = TRUE)

# create variable to store copy number
list_gr <- GRangesList()

# iterate over each patient, extract purity, ploidy and copy
# number segments
for(patient_id in patient_ids) {
  # set patient variables
  print(patient_id)
  patient_row <- which(patient_info$Patient == patient_id)
  path_patient <- file.path(tracerx_dir, patient_id)
  sex <- ifelse(unique(patient_info$Sex[patient_row]) == "Female", "XX", "XY")
  setwd(path_patient)

  # Get associated sample IDs
  sample_ids <- patient_info$Sample[patient_row]

    # iterate over sample IDs
    for(sample_id in sample_ids){
      # get sample row id
      print(sample_id)
      rw_id <-
        patient_info[, which(grepl(patient_id, Patient) & grepl(sample_id, Sample))]

      # Get single sample RRBS purity and ploidy 
      file_nm <- paste0(patient_id, ".", sample_id, ".ascat.output.RData")
      file_pth <-
        file.path(path_patient, "Copy_number", sample_id, "raw", file_nm)
      rho <- NA
      psi_grid <- NA
      psi_plot <- NA
      frac_LOH <- NA
      WGD <-  as.logical(NA)

      if (file.exists(file_pth)) {
        # load single sample ascat copy number segments
        load(file_pth)
        rm(file_pth)

        # if there was a unique solution, store raw single
        # sample resutls
        if(!is.null(ascat.output$segments)){
          rho <- as.numeric(ascat.output$aberrantcellfraction)
          psi_grid  <- as.numeric(ascat.output$ploidy)

          seg <- data.table(ascat.output$segments[, 2:6])
          seg[, length := endpos - startpos]

          sum_LOH <- seg[nMinor == 0, sum(length)]
          sum_not_LOH <- seg[nMinor != 0, sum(length)]
          frac_LOH <- sum_LOH / (sum_LOH + sum_not_LOH)

          # use PCAWG equation to determine WGD status
          WGD <- ifelse(2.95 - (2 * frac_LOH) - psi_grid < 0, TRUE, FALSE)
          psi_plot <- seg[, sum((nMajor + nMinor) * length) / sum(length)]

          if (patient_id %in% c("LTX143", "LTX164")) {
            counter <- length(list_gr) + 1

            list_gr[[counter]] <- 
              makeGRangesFromDataFrame(seg, 
                                      start.field = "startpos",
                                      end.field = "endpos",
                                      keep.extra.columns = TRUE)

            list_gr[[counter]]$WGD <-  WGD
            list_gr[[counter]]$hist <-  patient_info$Histology[rw_id]
            list_gr[[counter]]$patient <-  patient_id
            list_gr[[counter]]$sample <-  sample_id
          }
        }
      }

      # store single sample purity and ploidy values into data table
      patient_info$rho_raw[rw_id] <- rho
      patient_info$psi_grid_raw[rw_id] <- psi_grid
      patient_info$psi_plot_raw[rw_id] <- psi_plot
      patient_info$frac_LOH_raw[rw_id] <- frac_LOH
      patient_info$WGD_raw[rw_id] <- WGD

      # Get multi sample RRBS purity and ploidy
      if (!patient_id %in% c("LTX143", "LTX164")) {
        # load the resegmented copy number data
        file_nm <-
          paste0(patient_id, ".", sample_id,
                 ".ascat.output.resegmented.RData")
        file_pth <-
          file.path(path_patient, "Copy_number", sample_id, 
                    "re_segmented", file_nm)

        # set variables
        rho <- NA
        psi_grid <- NA
        psi_plot <- NA
        frac_LOH <- NA
        WGD <-  as.logical(NA)

        if (file.exists(file_pth)) {
          # load the multi-sample resegmented data
          load(file_pth)
          rm(file_pth)

          # extract the copy number segments, purity, ploidy and WGD status
          if (!is.null(ascat.output$segments)) {
            # extract the purity and ploidy
            counter <- length(list_gr) + 1
            rho <- as.numeric(ascat.output$aberrantcellfraction)
            psi_grid  <- as.numeric(ascat.output$ploidy)

            # get ASCN copy number segments
            seg <- data.table(ascat.output$segments[, 2:6])
            seg[,length := endpos-startpos]

            # convert to GRanges
            list_gr[[counter]] <-
              makeGRangesFromDataFrame(seg,
                                      start.field = "startpos",
                                      end.field = "endpos",
                                      keep.extra.columns = TRUE)

            sum_LOH <- seg[nMinor == 0, sum(length)]
            sum_not_LOH <- seg[nMinor != 0, sum(length)]
            frac_LOH <- sum_LOH / (sum_LOH + sum_not_LOH)

            WGD <- ifelse(2.95 - (2 * frac_LOH) - psi_grid < 0, TRUE, FALSE)

            list_gr[[counter]]$WGD <-  WGD
            list_gr[[counter]]$hist <- patient_info$Histology[rw_id]
            list_gr[[counter]]$patient <-  patient_id
            list_gr[[counter]]$sample <-  sample_id

            psi_plot <-
              seg[, sum((nMajor + nMinor) * length) / sum(length)]
          }
        }
      }
      # store multi-sample purity and ploidy values
      patient_info$rho_reseg[rw_id] <- rho
      patient_info$psi_grid_reseg[rw_id] <- psi_grid
      patient_info$psi_plot_reseg[rw_id] <- psi_plot
      patient_info$frac_LOH_reseg[rw_id] <- frac_LOH
      patient_info$WGD_reseg[rw_id] <- WGD

      # set WGS purity and ploidy variables
      file_nm <- paste0(sample_id, ".ascat.output.RData")
      file_pth <- file.path(path_patient, "WGS", sample_id, file_nm)
      rho <- NA
      psi_grid <- NA
      psi_plot <- NA
      frac_LOH <- NA
      WGD <-  as.logical(NA)

      if (file.exists(file_pth)) {
        # load WGS ascat copy number segments
        load(file_pth)
        rm(file_pth)
        if (!is.null(ascat.output$segments)) {
          rho <- as.numeric(ascat.output$aberrantcellfraction)
          psi_grid  <- as.numeric(ascat.output$ploidy)

          seg <- data.table(ascat.output$segments[, 2:6])
          seg[,length := endpos - startpos]

          sum_LOH <- seg[nMinor == 0, sum(length)]
          sum_not_LOH <- seg[nMinor != 0, sum(length)]
          frac_LOH <- sum_LOH / (sum_LOH + sum_not_LOH)

          WGD <-
            ifelse(2.95 - (2 * frac_LOH) - psi_grid < 0, TRUE, FALSE)

          psi_plot <-
            seg[, sum((nMajor + nMinor) * length) / sum(length)]
        }
      }

      # store WGS purity and ploidy values
      patient_info$WGS_rho[rw_id] <- rho
      patient_info$WGS_psi_grid[rw_id] <- psi_grid
      patient_info$WGS_psi_plot[rw_id] <- psi_plot
      patient_info$WGS_frac_LOH[rw_id] <- frac_LOH
      patient_info$WGS_WGD[rw_id] <- as.logical(WGD)

      # clean-up
      rm(
        seg, ascat.output, sum_LOH, sum_not_LOH,
        rho, psi_plot, psi_grid, frac_LOH, WGD
        )
   }
}

# reset work dir
setwd(work_dir)

# manual column re-ordering
id_cols <- c("Patient", "CRUKID", "Sample", "ID", "LTX_ID", "Histology", "Sex")
exome_cols <- c("Exome_psi", "Exome_rho", "VAF_rho")
rrbs_cols <-
  c("rho_raw", "psi_grid_raw", "psi_plot_raw", "frac_LOH_raw", "WGD_raw",
    "rho_reseg", "psi_grid_reseg", "psi_plot_reseg", "frac_LOH_reseg", "WGD_reseg")
wgs_cols <- c("WGS_rho", "WGS_psi_grid", "WGS_psi_plot",
              "WGS_frac_LOH", "WGS_WGD")
patient_info <-
  patient_info[, .SD, .SDcols = c(id_cols, exome_cols, rrbs_cols, wgs_cols)]
#"Published?", "immune_cluster", "cov_t", "cov_b", "ntrpcc", "CAMDAC_score"

################################################################################
# Annotate clinical info onto data table
################################################################################

# add additional clinical info
tmp <- fread(file.path(tracerx_dir, "tx100_RRBS_patients.csv"))
rws_match <- match(patient_info$Patient, tmp$SampleID)
patient_info$ethnicity <- tmp[rws_match, Ethnicity]
patient_info$age <- tmp[rws_match, Age]
patient_info$smoking <- tmp[rws_match, Smoking]
patient_info$pckyrs <- tmp[rws_match, Pyrs]
patient_info$OS <- tmp[rws_match, OS]
patient_info$RFS <- tmp[rws_match, RFS]
patient_info$Recurrence <- tmp[rws_match, Recurrence]
patient_info$Death <- tmp[rws_match, Death]
patient_info$Death_of_lung_cancer <- tmp[rws_match, Death_of_lung_cancer]
rm(tmp)

# removed failed run and lymph node samples
patient_info <-
  patient_info[paste(CRUKID, Sample, sep = "_") != "CRUK0084_R4", ]
patient_info <- patient_info[!grepl("LN", Sample),]

# print per sample and per patient WGD status
table(
  patient_info$WGD_reseg[patient_info$Histology == "Invasive adenocarcinoma"]
  )
table(
  patient_info$WGD_reseg[patient_info$Histology != "Invasive adenocarcinoma"]
  )
tmp <- patient_info[, .(mean(WGD_reseg)), by = .(Patient, Histology)]
table(tmp$V1[tmp$Histology == "Invasive adenocarcinoma"] > 0.5)
table(tmp$V1[tmp$Histology != "Invasive adenocarcinoma"] > 0.5)

################################################################################
# Annotate the normalised fraction LOH per patient and by AS-copy number
################################################################################

# custom overlaps
overlaps <- function(v, w, x, y) {
  value <- FALSE
  if (sum(!is.na(c(v, w, x, y))) == 4) {
    if (x <= v & v <= y) value <- TRUE
    if (x <= w & w <= y) value <- TRUE
  }
  return(value)
}

# combine GRanges
gr <- sort(unlist(list_gr))
gr$ID <- paste0(gr$patient, "_", gr$sample)
rm(list_gr)

# create object to store segments
counter <- 0
list_dfs <- list()
list_dfs_segment <- list()

# check if the file name exists
out_f_nm <- file.path(work_dir, in_dir, "frac_losses_and_gains.RData")

if (!file.exists(out_f_nm)) {
  for (ch in as.character(c(1:22))){
    counter <- counter + 1

    # get chrom gr
    gr_sub <- gr[as.character(seqnames(gr)) == ch]
    gr_sub$ASCN <- paste(gr_sub$nMajor, gr_sub$nMinor, sep = "_")

    # get min and max copy number segment coordinates
    minn <- min(start(gr_sub))
    maxx <- max(end(gr_sub))

    # get bin min to max sequence lenght out
    l <- length(seq(minn, maxx, 1E7))

    # create bins summary df
    df <- data.frame(chrom = ch,
                    start = seq(minn, maxx, 1E7)[1:(l - 1)],
                    end = seq(minn, maxx, 1E7)[2:l],
                    CN = numeric(length = l - 1),
                    frac_loss = numeric(length = l - 1),
                    frac_gain = numeric(length = l - 1),
                    frac_loss_LUAD = numeric(length = l - 1),
                    frac_gain_LUAD = numeric(length = l - 1),
                    frac_loss_LUSC = numeric(length = l - 1),
                    frac_gain_LUSC = numeric(length = l - 1),
                    frac_loss_LUAD_norm = numeric(length = l - 1),
                    frac_gain_LUAD_norm = numeric(length = l - 1),
                    frac_loss_LUSC_norm = numeric(length = l - 1),
                    frac_gain_LUSC_norm = numeric(length = l - 1),
                    stringsAsFactors = FALSE)

    # create bins df with per sample info
    df_segments <-
      data.frame(chrom = ch,
                start = seq(minn, maxx, 1E7)[1:(l - 1)],
                end = seq(minn, maxx, 1E7)[2:l],
                stringsAsFactors = FALSE)
    df_segments[, (1 + 3):(122 + 3)] <- as.numeric(NA)
    colnames(df_segments)[(1 + 3):(122 + 3)] <-
      as.character(patient_info$LTX_ID)
    df_segments_CN <- df_segments
    df_segments_nA <- df_segments
    rm(df_segments)

    # assign CN to each bin
    l <- length(gr_sub)
    for(i in 1:nrow(df)){
      flag <- mapply(function(w, x, y, z) overlaps(w, x, y, z),
                    w = rep(df$start[i], l),
                    x = rep(df$end[i], l),
                    y = start(gr_sub), z = end(gr_sub))
      gr_sub_sub <- gr_sub[flag]

      for(sampl in colnames(df_segments_CN)[(1+3):(122+3)]){
        gr_sub_sub_sub <-  gr_sub_sub[gr_sub_sub$ID==sampl]
        # get bin mean total, minor and major alleles copy number
        df_segments_CN[i, which(colnames(df_segments_CN) == sampl)] <-
          mean(gr_sub_sub_sub$nMajor + gr_sub_sub_sub$nMinor)
        df_segments_nA[i, which(colnames(df_segments_nA) == sampl)] <-
          mean(gr_sub_sub_sub$nMajor)
      }

      # get fractions of losses/gains
      df[i, "frac_loss"] <-
        sum(# LOH all involve losses
          grepl("_0", gr_sub_sub$ASCN) |
            # 2+1 is a loss in 4n genome
            (grepl("2_1", gr_sub_sub$ASCN) & gr_sub_sub$WGD == TRUE) |
            # 3+1 is a gain and a loss in 4n genome
            (grepl("3_1", gr_sub_sub$ASCN) & gr_sub_sub$WGD == TRUE)) /
              length(gr_sub_sub)
      df[i, "frac_gain"] <-
        sum(# 2+1 is a gain in 2n genome
          (grepl("2_1", gr_sub_sub$ASCN) & gr_sub_sub$WGD == FALSE) |
            # 3+1 is a gain for 2n and 4n genomes
            (grepl("3_1", gr_sub_sub$ASCN)) |
            # 2+0 is a loss and a gain in 2n genome
            (grepl("2_0", gr_sub_sub$ASCN) & gr_sub_sub$WGD == FALSE) |
            # 3+0 and above are gains and losses for 2n and 4n genomes
            (grepl("_0", gr_sub_sub$ASCN) &
              !grepl("1_0", gr_sub_sub$ASCN) & !grepl("2_0", gr_sub_sub$ASCN)) |
            # anything above copy numner 4 involved gains
            (gr_sub_sub$nMajor+gr_sub_sub$nMinor > 4)) /
              length(gr_sub_sub)

      gr_sub_sub_LUAD <- gr_sub_sub[grepl("adeno", gr_sub_sub$hist)]
      gr_sub_sub_LUSC <- gr_sub_sub[grepl("Squam", gr_sub_sub$hist)]

      df[i, "frac_loss_LUAD"] <-
        sum(# LOH all involve losses
          grepl("_0", gr_sub_sub_LUAD$ASCN) |
            # 2+1 is a loss in 4n genome
            (grepl("2_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == TRUE) |
            # 3+1 is a gain and a loss in 4n genome
            (grepl("3_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == TRUE)) /
              length(gr_sub_sub_LUAD)
      df[i, "frac_gain_LUAD"] <-
        sum(# 2+1 is a gain in 2n genome
          (grepl("2_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == FALSE) |
            # 3+1 is a gain for 2n and 4n genomes
            (grepl("3_1", gr_sub_sub_LUAD$ASCN)) |
            # 2+0 is a loss and a gain in 2n genome
            (grepl("2_0", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == FALSE) |
            # 3+0 and above are gains and losses for 2n and 4n genomes
            (grepl("_0", gr_sub_sub_LUAD$ASCN) &
              !grepl("1_0", gr_sub_sub_LUAD$ASCN) &
                !grepl("2_0", gr_sub_sub_LUAD$ASCN)) |
            # anything above copy numner 4 involved gains
            (gr_sub_sub_LUAD$nMajor+gr_sub_sub_LUAD$nMinor > 4)) /
              length(gr_sub_sub_LUAD)

      df[i, "frac_loss_LUSC"] <-
        sum(# LOH all involve losses
            grepl("_0", gr_sub_sub_LUSC$ASCN) |
            # 2+1 is a loss in 4n genome
            (grepl("2_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == TRUE) |
            # 3+1 is a gain and a loss in 4n genome
            (grepl("3_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == TRUE)) /
              length(gr_sub_sub_LUSC)
      df[i, "frac_gain_LUSC"] <-
        sum(# 2+1 is a gain in 2n genome
            (grepl("2_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == FALSE) |
            # 3+1 is a gain for 2n and 4n genomes
            (grepl("3_1", gr_sub_sub_LUSC$ASCN)) |
            # 2+0 is a loss and a gain in 2n genome
            (grepl("2_0", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == FALSE) |
            # 3+0 and above are gains and losses for 2n and 4n genomes
            (grepl("_0", gr_sub_sub_LUSC$ASCN) &
              !grepl("1_0", gr_sub_sub_LUSC$ASCN) &
                !grepl("2_0", gr_sub_sub_LUSC$ASCN)) |
            # anything above copy numner 4 involved gains
            (gr_sub_sub_LUSC$nMajor+gr_sub_sub_LUSC$nMinor > 4)) /
              length(gr_sub_sub_LUSC)

      tmp <- data.table(data.frame(gr_sub_sub))
      tmp <- tmp[, .(start = mean(start), end = mean(end), 
                    nMajor = round2(mean(nMajor), digits=0),
                    nMinor = round2(mean(nMinor), digits=0),
                    WGD = as.logical(round2(mean(as.numeric(WGD)), digits = 0))),
                by = .(seqnames, patient, hist)]
      gr_sub_sub <- makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE)
      gr_sub_sub$ASCN <- paste(gr_sub_sub$nMajor, gr_sub_sub$nMinor, sep = "_")
      rm(tmp)
      gr_sub_sub_LUAD <- gr_sub_sub[grepl("adeno", gr_sub_sub$hist)]
      gr_sub_sub_LUSC <- gr_sub_sub[grepl("Squam", gr_sub_sub$hist)]

      df[i, "frac_loss_LUAD_norm"] <-
        sum(# LOH all involve losses
          grepl("_0", gr_sub_sub_LUAD$ASCN) |
            # 2+1 is a loss in 4n genome
            (grepl("2_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == TRUE) |
            # 3+1 is a gain and a loss in 4n genome
            (grepl("3_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == TRUE)) /
              length(gr_sub_sub_LUAD)
      df[i, "frac_gain_LUAD_norm"] <-
        sum(# 2+1 is a gain in 2n genome
          (grepl("2_1", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == FALSE) |
            # 3+1 is a gain for 2n and 4n genomes
            (grepl("3_1", gr_sub_sub_LUAD$ASCN)) |
            # 2+0 is a loss and a gain in 2n genome
            (grepl("2_0", gr_sub_sub_LUAD$ASCN) & gr_sub_sub_LUAD$WGD == FALSE) |
            # 3+0 and above are gains and losses for 2n and 4n genomes
            (grepl("_0", gr_sub_sub_LUAD$ASCN) &
              !grepl("1_0", gr_sub_sub_LUAD$ASCN) &
                !grepl("2_0", gr_sub_sub_LUAD$ASCN)) |
            # anything above copy numner 4 involved gains
            (gr_sub_sub_LUAD$nMajor+gr_sub_sub_LUAD$nMinor > 4)) /
              length(gr_sub_sub_LUAD)

      df[i, "frac_loss_LUSC_norm"] <-
        sum(# LOH all involve losses
          grepl("_0", gr_sub_sub_LUSC$ASCN) |
            # 2+1 is a loss in 4n genome
            (grepl("2_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == TRUE) |
            # 3+1 is a gain and a loss in 4n genome
            (grepl("3_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == TRUE)) /
              length(gr_sub_sub_LUSC)
      df[i, "frac_gain_LUSC_norm"] <-
        sum(# 2+1 is a gain in 2n genome
          (grepl("2_1", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == FALSE) |
            # 3+1 is a gain for 2n and 4n genomes
            (grepl("3_1", gr_sub_sub_LUSC$ASCN)) |
            # 2+0 is a loss and a gain in 2n genome
            (grepl("2_0", gr_sub_sub_LUSC$ASCN) & gr_sub_sub_LUSC$WGD == FALSE) |
            # 3+0 and above are gains and losses for 2n and 4n genomes
            (grepl("_0", gr_sub_sub_LUSC$ASCN) &
              !grepl("1_0", gr_sub_sub_LUSC$ASCN) &
                !grepl("2_0", gr_sub_sub_LUSC$ASCN)) |
            # anything above copy numner 4 involved gains
            (gr_sub_sub_LUSC$nMajor+gr_sub_sub_LUSC$nMinor > 4)) /
              length(gr_sub_sub_LUSC)
    }
    list_dfs[[counter]] <- df
    list_dfs_segment[[counter]] <- rbind(df_segments_CN, df_segments_nA)
  }
  # unlist results
  df <- bind_rows(list_dfs)
  df_segments <- bind_rows(list_dfs_segment)

  # save data
  save(df, df_segments, file = out_f_nm)

} else {
  load(out_f_nm)
}

# get chromosome start and end location for plot
df_plot <-
  data.table(df)[, .(start = min(which(df$chrom == chrom)),
                     end = max(which(df$chrom == chrom)),
                     mid = max(which(df$chrom == chrom)) -
                           (max(which(df$chrom == chrom)) -
                            min(which(df$chrom == chrom))) / 2),
                  by = .(chrom)]

# plot fraction gains and losses LUAD
p_LUAD <- ggplot(df) + theme_minimal() + ylim(c(-1, 1.15)) + 
  geom_col(aes(x=1:nrow(df), y=-frac_loss_LUAD_norm), col="lightblue", fill="lightblue") +
  geom_col(aes(x=1:nrow(df), y=frac_gain_LUAD_norm), col="red3", fill="red3") +
  geom_text(data=df_plot, aes(x=mid, y=1.1, label=chrom), size=2) +
  geom_vline(data=df_plot, aes(xintercept=end+0.5), col="grey70", size=0.5)

ggsave(p_LUAD, 
       filename = file.path(work_dir, out_dir,
       "gains_and_losses_summary_LUAD_normalised.pdf"),
       device = "pdf", units = "in", width = 5.5, height = 2)

# plot fraction gains and losses LUSC
p_LUSC <- ggplot(df) + theme_minimal() + ylim(c(-1, 1.15))+ 
  geom_col(aes(x=1:nrow(df), y=-frac_loss_LUSC_norm), col="lightblue", fill="lightblue") +
  geom_col(aes(x=1:nrow(df), y=frac_gain_LUSC_norm), col="red3", fill="red3") +
  geom_text(data=df_plot, aes(x=mid, y=1.1, label=chrom), size=2) +
  geom_vline(data=df_plot, aes(xintercept=end+0.5), col="grey70", size=0.5)

ggsave(p_LUSC,
       filename = file.path(work_dir, out_dir,
       "gains_and_losses_summary_LUSC_normalised.pdf"),
       device = "pdf", units = "in", width = 5.5, height = 2)

################################################################################
# Compare ASCAT(.m) purity and ploidy across platforms
################################################################################

# load WGBS data
dt_prc <-
  fread(file.path(work_dir, in_dir, "mcrpc_wgbs_wgs_overlap.csv"),
        sep = ",",
        colClasses = c(perc_overlap = "numeric"))

# remove DTB-255 as the WGS solution as a large homozygous deletion and
# we do not have access to the raw data to refit this samples
dt_prc <- dt_prc[sample_id != "DTB-255", ]

# load WGBS data
dt_ccrc_rrbs <-
  fread(file.path(work_dir, in_dir, "CRC_RRBS_WGS_overlap.txt"))

# load WGBS data
dt_ccrc_wgbs <-
  fread(file.path(work_dir, in_dir, "CRC_WGBS_WGS_overlap.txt"))

# combine all copy number data
dt_plot <- rbind(
  patient_info[, .(dnaseq_rho = Exome_rho, bsseq_rho = rho_reseg, 
                   dnaseq_psi = Exome_psi, bsseq_psi = psi_grid_reseg,
                   method = "RRBS", shaps = "WES")],
  patient_info[!is.na(WGS_rho),
               .(dnaseq_rho = WGS_rho, bsseq_rho = rho_reseg,
                 dnaseq_psi = WGS_psi_grid, bsseq_psi = psi_grid_reseg,
                 method = "RRBS", shaps = "WGS")],
  dt_prc[, .(dnaseq_rho = purity_a[tolower(name_a) == "wgs"],
             bsseq_rho = purity_b[tolower(name_b) == "wgbs"],
             dnaseq_psi = ploidy_a[tolower(name_a) == "wgs"],
             bsseq_psi = ploidy_b[tolower(name_b) == "wgbs"],
             method = "WGBS", shaps = "WGS")],
  dt_ccrc_rrbs[, .(dnaseq_rho = purity_a[tolower(name_a) == "wgs"],
              bsseq_rho = purity_b[tolower(name_b) == "rrbs"],
              dnaseq_psi = ploidy_a[tolower(name_a) == "wgs"],
              bsseq_psi = ploidy_b[tolower(name_b) == "rrbs"],
              method = "RRBS", shaps = "WGS")],
  dt_ccrc_wgbs[, .(dnaseq_rho = purity_a[tolower(name_a) == "wgs"],
              bsseq_rho = purity_b[tolower(name_b) == "wgbs"],
              dnaseq_psi = ploidy_a[tolower(name_a) == "wgs"],
              bsseq_psi = ploidy_b[tolower(name_b) == "wgbs"],
              method = "WGBS", shaps = "WGS")]
)
dt_plot[, comparison := paste(method, "vs", shaps)]
dt_plot[, table(comparison)]

# plot rho BS-seq vs other platforms
dt_plot[, cor(dnaseq_rho, bsseq_rho), by = .(shaps)]
out_f_nm <-
  file.path(work_dir, out_dir,
            "/bsseq_vs_dnaseq_purity_all_samples.pdf")
p <- ggplot(dt_plot, aes(x = dnaseq_rho, y = bsseq_rho)) +
  geom_point(aes(color = comparison, shape = comparison), size = 3,
             alpha = ifelse(dt_plot[, comparison == "RRBS vs WGS"], 1, 0.6)) +
  theme_classic() +
  scale_color_manual(name = "comparison",
                     values = c("RRBS vs WES" = "cornflowerblue",
                                "RRBS vs WGS" = "black",
                                "WGBS vs WGS" = "firebrick")) +
  scale_shape_manual(name = "comparison",
    values = c("RRBS vs WES" = 15, "RRBS vs WGS" = 16, "WGBS vs WGS" = 18)) +
  scale_y_continuous("Matched DNA-Seq purity estimate",
    breaks = seq(0.1, 0.9, 0.2), limits = c(.1, .9)) +
  scale_x_continuous("BS-Seq purity estimate",
    breaks = seq(0.1, 0.9, 0.2), limits = c(.1, .9))
ggsave(p, filename = out_f_nm,
       device = "pdf", units = "in", width = 4.5, height = 3)

# plot ploidy BS-seq vs other platforms
dt_plot[, cor(dnaseq_psi, bsseq_psi), by = .(shaps)]
p_psi <- ggplot(dt_plot, aes(x = dnaseq_psi, y = bsseq_psi)) +
  geom_point(aes(color = comparison, shape = comparison),
             size = 3, alpha = ifelse(dt_plot[, comparison == "RRBS vs WGS"], 1, 0.6)) +
  theme_classic() + 
  scale_color_manual(name="comparison",
                    values = c("RRBS vs WES" = "cornflowerblue",
                               "RRBS vs WGS" = "black",
                               "WGBS vs WGS" = "firebrick")) +
  scale_shape_manual(name="comparison",
    values = c("RRBS vs WES" = 15, "RRBS vs WGS" = 16, "WGBS vs WGS" = 18)) +
  scale_y_continuous("Matched DNA-Seq ploidy estimate",
    breaks = seq(2, 5, 1), limits = c(1.5, 5.5)) +
  scale_x_continuous("BS-Seq ploidy estimate",
    breaks = seq(2, 5, 1), limits = c(1.5, 5.5))

# save
out_f_nm <-
  file.path(work_dir, out_dir, "/bsseq_vs_dnaseq_ploidy_all_samples.pdf")
ggsave(p_psi, filename = out_f_nm,
       device = "pdf", units = "in", width = 4.5, height = 3)

# plot rho vs rho_reseg
p_multi_vs_single <-
  ggplot(patient_info, aes(x = rho_reseg, y = rho_raw)) +
  theme_minimal() +
  geom_point(col = "grey30", size = 3, shape = 20) +
  xlab("CAMDAC purity multi-sample") +
  ylab("CAMDAC purity single sample")

# save
ggsave(p_multi_vs_single,
       filename = file.path(work_dir, out_dir, "purity_all_samples.pdf"),
       device = "pdf", units = "in", width = 3, height = 3)

# plot rho_reseg LUAD vs LUSC
tmp <- patient_info
tmp$hist <- ifelse(grepl("adeno", tmp$Histology), "LUAD", "LUSC")
p <-
  ggplot(tmp, aes(x = hist, y = rho_reseg, col = hist)) +
  theme_minimal() +
  geom_violin(position = "dodge", alpha = 0)  +
  geom_jitter(aes(fill = hist), alpha = .5, width = .1) +
  xlab("Histology") +
  ylab("ASCAT.m purity")

# save
ggsave(p, filename = file.path(work_dir, out_dir, "purity_by_hist.pdf"),
       device = "pdf", units = "in", width = 6, height = 3)

# plot psi LUAD vs LUSC
p <- ggplot(tmp, aes(x = hist, y = psi_grid_reseg, col = hist)) +
  theme_minimal() +
  geom_violin(position = "dodge", alpha = 0) +
  geom_jitter(aes(fill = hist), alpha = .5, width = .1) +
  xlab("Histology") +
  ylab("ASCAT.m ploidy")

# save
ggsave(p,
       filename = file.path(work_dir, out_dir, "ploidy_by_hist.pdf"),
       device = "pdf", units = "in", width = 4, height = 3)

# extract multi-sample ASCAT.m ploidy estimates, fraction of the
# genome with LOH and status
tmp <- data.frame(psi_grid = patient_info$psi_grid_reseg,
                  frac_LOH = patient_info$frac_LOH_reseg,
                  WGD = factor(as.logical(patient_info$WGD_reseg),
                               levels = c("TRUE", "FALSE"), ordered = TRUE),
                  patient = patient_info$Patient,
                  stringsAsFactors = FALSE)
tmp$histology <-
  as.character(patient_info$Histology[match(tmp$patient, patient_info$Patient)])
tmp$histology <-
  ifelse(as.character(tmp$histology) == "Invasive adenocarcinoma",
         "LUAD", "LUSC")
tmp$shapes <-
  ifelse(tmp$histology == "LUAD", 16,
    ifelse(tmp$histology == "LUSC", 18, 17))

# plot psi_grid vs frac_LOH
p_frac_LOH <- ggplot(tmp) +
  geom_point(aes(x = frac_LOH, y = psi_grid, col = WGD), size = 3, alpha = 0.9,
             shape = tmp$shapes) +
  xlab("ASCAT.m fraction LOH") +
  theme_classic() +
  scale_y_continuous(name = "ASCAT.m ploidy",
    limits = c(1,5.25), breaks = seq(1,5,1)) +
  scale_color_manual(values = c("FALSE" = "gold2", "TRUE" = "deeppink3")) +
  geom_abline(slope = -2, intercept = 2.95, col = "grey50",
    linetype = "dashed", size = .5)

# save
out_f_nm <- file.path(work_dir, out_dir, "frac_LOH_vs_psi_grid_all_samples.pdf")
ggsave(p_frac_LOH, filename = out_f_nm, 
       device = "pdf", units = "in", width = 3.75, height = 3)

################################################################################
# Compare the allele specific copy number overlap across platforms
################################################################################

# Read in tracerx data
load("inputs/dt_tracerx.RData")
dt_tracerx <- as.data.table(dt_tracerx)
dt_tracerx[, "perc_overlap" := as.numeric(as.character(perc_copy_overlap))]
dt_ccrc_rrbs[, "perc_overlap" := as.numeric(as.character(perc_overlap))]
dt_ccrc_wgbs[, "perc_overlap" := as.numeric(as.character(perc_overlap))]

# combine all copy number data
dt_with_matched_normal <- rbind(
  dt_tracerx[, .(perc_overlap, "perc_overlap_1" = perc_copy_within_1,
    method = "RRBS", shaps = "WES")],
  dt_ccrc_rrbs[, .(perc_overlap, perc_overlap_1, method = "RRBS", shaps = "WGS")],
  dt_ccrc_wgbs[, .(perc_overlap, perc_overlap_1, method = "WGBS", shaps = "WGS")]          
)
dt_with_matched_normal[, comparison := paste(method, "vs", shaps)]
dt_with_matched_normal[, table(comparison)]

# plot ploidy BS-seq vs other platforms
dt_with_matched_normal[, "x" :=
  nrow(dt_with_matched_normal) - rank(perc_overlap)
  ]
p_matched_normal <-
  ggplot(dt_with_matched_normal) +
  geom_point(aes(x = x, y = perc_overlap, shape = comparison),
             size = 3, alpha = 0.95, color = "orange2") +
  geom_point(aes(x = x, y = perc_overlap_1, shape = comparison),
             size = 3, alpha = 0.95, color = "darkgreen") +
  theme_classic() +
  scale_shape_manual(name = "comparison",
    values = c("RRBS vs WES" = 15, "RRBS vs WGS" = 16, "WGBS vs WGS" = 18)) +
  scale_y_continuous("Percentage overlap", breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_x_continuous("Sample", breaks = seq(0, nrow(dt_with_matched_normal), 25)) 

# save
out_f_nm <- file.path(work_dir, out_dir, "bsseq_vs_dnaseq_ascn_overlap.pdf")
ggsave(p_matched_normal, filename = out_f_nm,
       device = "pdf", units = "in", width = 5, height = 3)

# plot % overlap
dt_prc[, "x" := nrow(dt_prc) - rank(perc_overlap)]
dt_prc[, "comparison" := "WGBS vs WGS"]
p_unmatched_normal <-
  ggplot(dt_prc) +
  geom_point(aes(x = x, y = perc_overlap, shape = comparison),
             size = 3, alpha = 0.9, color = "orange2") +
  geom_point(aes(x = x, y = perc_overlap_1, shape = comparison),
             size = 3, alpha = 0.9, color = "darkgreen") +
  theme_classic() +
  scale_shape_manual(name = "comparison", values = c("WGBS vs WGS" = 18)) +
  scale_y_continuous("Percentage overlap", breaks = seq(0, 100, 10), limits = c(0, 100)) +
  scale_x_continuous("Sample", breaks = seq(0, 80, 10), limits = c(0, 80))

# save
out_f_nm <-
  file.path(work_dir, out_dir,
            "bsseq_vs_dnaseq_ascn_overlap_panel_of_normals.pdf")
ggsave(p_unmatched_normal, filename = out_f_nm,
       device = "pdf", units = "in", width = 4.5, height = 3)

# END