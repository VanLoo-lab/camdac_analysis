#' @description Get sample pairwise distance
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
x <- c("stringr", "grid", "gridExtra", "ggdendro", "scales", "ggplot2",
       "data.table", "readr", "ggfortify", "ggridges")
load_lib <- function (x) {
  invisible(
    lapply(x, function(y) {
      suppressWarnings(
        suppressMessages(
          library(y, character.only = TRUE, quietly = TRUE,
                  warn.conflicts = FALSE)
          )
        )
      }
    )
  )
}
load_lib(x)
rm(x)

################################################################################
# Source useful functions
################################################################################

#' @description Get sample pairwise distance
#' @param dt a data.table
pearson <- function(dt, ...) {
  nc <- ncol(dt)

  tmp <- dt[, list(list(1 - cor(dt, method = "pearson", use = "everything")))]
  mat <- matrix(unlist(tmp), nrow = nc)

  rownames(mat) <-  colnames(dt)
  colnames(mat) <-  colnames(dt)

  mat <- as.dist(mat)
  attr(mat, "distfun") <- "pearson"
  attr(mat, "Labels") <- colnames(dt)

  return(mat)
}

################################################################################
# Load patient information
################################################################################

# set paths
parent_dir <- "~/rdeocd/txMeth"
parent_dir_facs <-  "~/rdeocd/FACS_sorted_txMeth"
work_dir <- "~/data/CAMDAC_paper/figure2N"

# load FACS patient info
f_nm <- file.path(parent_dir_facs, "data",
  "NuGEN_summary_report_FACS_sorted_samples.txt")
patient_info_facs <- fread(file = f_nm, sep = "\t")

# annotate the population (1 = normal, 2 = tumor)
patient_info_facs[(nchar(as.character(Patient)) == 6L & substr(as.character(Patient), 4, 4) == "0"),
  pop_id := substr(Sample_Name, 7, 10)]
patient_info_facs[!(nchar(as.character(Patient)) == 6L & substr(as.character(Patient), 4, 4) == "0"),
  pop_id := substr(Sample_Name, 8, 11)]

# remove any "," in the numbers
patient_info_facs <-
  patient_info_facs[, lapply(.SD,
    function(x) gsub(pattern = ",", replacement = "", x)
    ), .SDcols = colnames(patient_info_facs)]

# convert the column to numeric
col_num <- colnames(patient_info_facs)[7:(ncol(patient_info_facs) - 1)]
patient_info_facs[, as.character(col_num) :=
  lapply(.SD, as.numeric), .SDcols = col_num]

# load FACS patient info
f_nm <-
  paste0(parent_dir,
         "/tx100clinicalinfo_RRBS_published_status_",
         "immune_clusters_purity_ploidy_CpG_cov.txt")
patient_info <- fread(f_nm, sep = "\t")
rm(f_nm)

# remove lymph nodes
patient_info <- patient_info[!grepl("LN", Sample), ]

# extract the patient IDs of all TRACERx patients
patient_ids <- patient_info[, unique(Patient)]

# annotate CRUKID onto FACS data.table
crukids <- patient_info[match(patient_info_facs$Patient, Patient), CRUKID]
patient_info_facs[, CRUKID := crukids]

# annotate purity and ploidy onto FACS data.table
ltxids_facs <- patient_info_facs[, paste0(Patient, "_", Region)]
rho_vals <- patient_info[match(ltxids_facs, LTX_ID), rho_reseg]
patient_info_facs[, rho := rho_vals]
psi_vals <- patient_info[match(ltxids_facs, LTX_ID), psi_reseg]
patient_info_facs[, psi := psi_vals]

################################################################################
# Load bulk and CAMDAC-purified tumours for the TRACERx lung data
################################################################################

# extract the IDs of the patients with FACS data
patient_ids <- patient_info[, unique(Patient)]

# create merge function
my_merge <-
  function(...) merge(..., by = c("chromosome", "start", "end"), all = FALSE)

# create patient specific results object for merging
list_patientss_mt <- list()
pos_patientss_mt <- NULL
list_patientss_mb <- list()
pos_patientss_mb <- NULL

# get data for each patient
for (patient_id in patient_ids) {
  # set path
  path_patient <- file.path(parent_dir, patient_id)
  setwd(path_patient)

  # Create object to store per patient output
  list_patient_mt <- list()
  list_patient_mb <- list()
  cols <- c("chromosome", "start", "end", "m_b", "m_t_corr", "m_n")
  # each list element contains the data for a given sample of a given patient
  # The data in each list element is that described by cols

  # Get sample ID and CRUK ID for the ith patient
  sample_ids <- as.character(
    patient_info[which(Patient == patient_id), Sample]
  )

  # get data for each sample
  for (sample_id in sample_ids){
    i <- which(sample_ids == sample_id)
    cat(paste(patient_id, sample_id, sep = "_"))

    # check deconvoluted methylations rate estimates files exists and load file
    f_nm <- file.path(path_patient, "Methylation", sample_id, "re_segmented",
      "purified_tumour_DMP_calls.RData")
    if (!file.exists(f_nm)) {
      f_nm <- file.path(path_patient, "Methylation", sample_id, "raw",
        "purified_tumour_DMP_calls.RData")
        if (!file.exists(f_nm)) {
          f_nm <- file.path(path_patient, "Methylation", sample_id,
            "purified_tumour_DMP_calls.RData")
        }
      load(f_nm)
    } else if (file.exists(f_nm)) {
      load(f_nm)
    }
    rm(f_nm)

    # load methylation info
    cat(" - bulk and deconvoluted methylation rates loaded", "\n")
    df_mt <- data.table(df_mt)
    df_mt[, "chromosome" := as.character(chromosome)]

    # get label
    label <- paste0(patient_id, "_", sample_id)

    # store methylation info and CpG loci
    if (i != length(sample_ids)) {
      # get mt
      tmp <- df_mt
      list_patient_mt[[length(list_patient_mt) + 1]] <-
        tmp[, .SD, .SDcols = cols[c(1:3, 5)]]
      colnames(list_patient_mt[[length(list_patient_mt)]])[4] <- label

      # get mb
      list_patient_mb[[length(list_patient_mb) + 1]] <-
        tmp[, .SD, .SDcols = cols[c(1:4)]]
      colnames(list_patient_mb[[length(list_patient_mb)]])[4] <- label
    }

    # store normal methylation info
    if (i == length(sample_ids)) {
      # get mt
      tmp <- df_mt
      list_patient_mt[[length(list_patient_mt) + 1]] <-
        tmp[, .SD, .SDcols = cols[c(1:3, 5:6)]]
      colnames(list_patient_mt[[length(list_patient_mt)]])[4:5] <-
        c(label, paste0(patient_id, "_Normal"))

      # get mb
      list_patient_mb[[length(list_patient_mb) + 1]] <-
        tmp[, .SD, .SDcols = cols[c(1:4, 6)]]
      colnames(list_patient_mb[[length(list_patient_mb)]])[4:5] <-
        c(label, paste0(patient_id, "_Normal"))
    }

    # get pos DMP mt
    tmp <-
      df_mt[abs(m_t_corr - m_n) > 0.2 &
            (prob_hypo < 0.01 | prob_hyper < 0.01), ]
    pos_patient_mt <- tmp[, .SD, .SDcols = cols[c(1:3)]]

    # get pos DMP mb
    tmp <-
      df_mt[abs(m_b - m_n) > 0.2 &
            (prob_hypo < 0.01 | prob_hyper < 0.01), ]
    pos_patient_mb <- tmp[, .SD, .SDcols = cols[c(1:3)]]

    rm(df_mt)
  }

  list_patientss_mt[[length(list_patientss_mt) + 1]] <-
    Reduce(my_merge, list_patient_mt)
  list_patientss_mb[[length(list_patientss_mb) + 1]] <-
    Reduce(my_merge, list_patient_mb)

  if (!is.null(pos_patientss_mt)) {
    tmpt <- rbindlist(list(pos_patient_mt, pos_patientss_mt))
    tmpb <- rbindlist(list(pos_patient_mb, pos_patientss_mb))
  } else {
    tmpb <- pos_patient_mb
    tmpt <- pos_patient_mt
  }

  # remove duplicate entries
  tmpb <- tmpb[!duplicated(tmpb)]
  tmpt <- tmpt[!duplicated(tmpt)]

  pos_patientss_mt <- tmpt
  pos_patientss_mb <- tmpb

  rm(list_patient_mt, list_patient_mb, i, sample_ids,
     pos_patient_mt, pos_patient_mb, tmp, tmpb, tmpt)
}

################################################################################
# Load FACS-purified tumours
################################################################################

# create patient specific results object for merging
list_patientss_mf <- list()

for (sample_id in patient_info_facs[, sort(Sample_Name)]) {
  patient_id <- patient_info_facs[Sample_Name == sample_id, Patient]
  region_id <- patient_info_facs[Sample_Name == sample_id, Region]
  pop_id <- paste0("pop", patient_info_facs[Sample_Name == sample_id, Pop])
  print(sample_id)

  # create dirs
  path_patient <- file.path(parent_dir_facs, "data", patient_id)
  setwd(path_patient)

  cols <- c("chromosome", "start", "end", "m_t_corr", "m_n")
  # each list element contains the data for a given sample of a given patient
  # The data in each list element is that described by cols

  cat(paste(patient_id, sample_id, sep = "_"))

  # check deconvoluted methylations rate estimates files exists and load file
  f_nm <- file.path(path_patient, "Methylation", sample_id,
    "purified_tumour_DMP_calls.RData")
  load(f_nm)

  # load methylation info
  cat(" - facs deconvoluted methylation rates loaded", "\n")
  df_mt <- data.table(df_mt)
  df_mt[, "chromosome" := as.character(chromosome)]
  df_mt[, "m_t_corr" := m_b]

  # store methylation info and CpG loci
  if (pop_id == "pop2") {
    # get mt
    list_patientss_mf[[length(list_patientss_mf) + 1]] <-
      df_mt[, .SD, .SDcols = c(cols[c(1:4)])]

    # get label
    label <- paste0(patient_id, "_", region_id)
    colnames(list_patientss_mf[[length(list_patientss_mf)]])[4] <- label
  }

  # store normal methylation info
  if (pop_id == "pop1") {
    # get mn
    list_patientss_mf[[length(list_patientss_mf) + 1]] <-
        df_mt[, .SD, .SDcols = c(cols[c(1:3, 5)])]

    # get label
    label <- paste0(patient_id, "_", region_id, "_Normal")
    colnames(list_patientss_mf[[length(list_patientss_mf)]])[4] <- label
  }
}

################################################################################
# Compute the correlations
################################################################################

# concatenate methylation info at all CpG loci and removing NAs
dt_mb_all <- Reduce(my_merge, list_patientss_mb)
dt_mb_all <- dt_mb_all[rowSums(!is.na(dt_mb_all)) == ncol(dt_mb_all), ]
dt_mb_all[, ID := paste(chromosome, start, end, sep = "_")]
dt_mb_all <- dt_mb_all[!duplicated(ID), ]
rm(list_patientss_mb)

# concatenate methylation info at all CpG loci and removing NAs
dt_mt_all <- Reduce(my_merge, list_patientss_mt)
dt_mt_all <- dt_mt_all[rowSums(!is.na(dt_mt_all)) == ncol(dt_mt_all), ]
dt_mt_all[, ID := paste(chromosome, start, end, sep = "_")]
dt_mt_all <- dt_mt_all[!duplicated(ID), ]
rm(list_patientss_mt)

# concatenate methylation info at all CpG loci and removing NAs
dt_mf_all <- Reduce(my_merge, list_patientss_mf)
dt_mf_all <- dt_mf_all[rowSums(!is.na(dt_mf_all)) == ncol(dt_mf_all), ]
dt_mf_all[, ID := paste(chromosome, start, end, sep = "_")]
dt_mf_all <- dt_mf_all[!duplicated(ID), ]
rm(list_patientss_mf)

# get CpGs with no NAs in mb, mt or mf set but with >0 DMP from CAMDAC mt
pos_patientss_mt[, ID := paste(chromosome, start, end, sep = "_")]
pos_patientss_mt <- pos_patientss_mt[!duplicated(ID)]
pos_patientss_mt <- pos_patientss_mt[ID %in% dt_mf_all$ID, ]
pos_patientss_mt <- pos_patientss_mt[ID %in% dt_mt_all$ID, ]
pos_patientss_mt[, ID := NULL]

# reduce dt to CpG that are DMP loci in one or more sample
dt_mf <- Reduce(my_merge, list(dt_mf_all, pos_patientss_mt))
dt_mf <- dt_mf[, 4:ncol(dt_mf)]
dt_mf[, ID := NULL]
rm(dt_mf_all)

# reduce dt to CpG that are DMP loci in one or more sample
dt_mt <- Reduce(my_merge, list(dt_mt_all, pos_patientss_mt))
dt_mt <- dt_mt[, 4:ncol(dt_mt)]
dt_mt[, ID := NULL]
rm(dt_mt_all)

# reduce dt to CpG that are DMP loci in one or more sample
dt_mb <- Reduce(my_merge, list(dt_mb_all, pos_patientss_mt))
dt_mb <- dt_mb[, 4:ncol(dt_mb)]
dt_mb[, ID := NULL]
rm(dt_mb_all)

# save data
setwd(work_dir)
saveRDS(dt_mb, file = "inputs/lung_mb.RDS")
saveRDS(dt_mt, file = "inputs/lung_mt.RDS")
saveRDS(dt_mf, file = "inputs/lung_mf.RDS")

# FACS normal vs. bulk normal and FACS-purified tumour
# vs. CAMDAC-purified tumour
cols <- colnames(dt_mf)
dist_mtf <-
  data.table(sample = cols,
             x = 0,
             flag = "flag",
             flag2 = "FACS")
for (coll in cols) {
  if (grepl("Normal", coll)) {
    coll2 <- paste0(substr(coll, 1, 6), "_Normal")
    flagg <- "nn"
  } else {
    coll2 <- coll
    flagg <- "tt_intra"
  }
  print(coll)
  print(coll2)
  y <- cor(x = dt_mf[, .SD, .SDcols = coll],
           y = dt_mt[, .SD, .SDcols = coll2], method = "pearson")
  y <- 1 - y
  dist_mtf[sample == coll, "x" := y]
  dist_mtf[sample == coll, "flag" := flagg]
  if (flagg == "tt_intra") {
    dist_mtf_tmp <-
      data.table(sample = coll,
             x = 0,
             flag = "flag",
             flag2 = "FACS")
    y <- cor(x = dt_mf[, .SD, .SDcols = coll],
             y = dt_mb[, .SD, .SDcols = coll2], method = "pearson")
    y <- 1 - y
    dist_mtf_tmp[sample == coll, "x" := y]
    dist_mtf_tmp[sample == coll, "flag" := "bb_intra"]
    dist_mtf <- rbind(dist_mtf, dist_mtf_tmp)
    rm(dist_mtf_tmp)
  }
}

# CAMDAC-purified tumour sample pair correlations
# vs. bulk sample pair correlations
for (name in c("mb", "mt")) {
  dt <- get(ls(pattern = paste("dt", name, sep = "_")))

  # scale data.table
  cols <- colnames(dt)
  dt[, (cols) := lapply(.SD, scale), .SDcols = cols]
  rm(cols)

  # make distance matrix
  dista <- pearson(dt)
  assign(paste0("dist_", name), dista)
  rm(dista)
}

# male vs female
mal <- as.character(patient_info[Sex == "Male", unique(Patient)])
fem <- as.character(patient_info[Sex == "Female", unique(Patient)])
fem <- fem[fem != resolve_tx_patient_id("CRUK0047", .tx_map)]

# focus on normal-normal
tmp <- as.data.frame(as.matrix(dist_mb))
tmp <- tmp[grepl("Normal", rownames(tmp)) & rownames(tmp) != paste0(resolve_tx_patient_id("CRUK0047", .tx_map), "_Normal"),
           grepl("Normal", colnames(tmp)) & colnames(tmp) != paste0(resolve_tx_patient_id("CRUK0047", .tx_map), "_Normal")]
dist_nn <- as.numeric(unlist(as.dist(as.matrix(tmp))))
tmp <- tmp[substr(rownames(tmp), 1, 6) %in% mal,
      substr(rownames(tmp), 1, 6) %in% fem]
dist_nn_fm <- as.numeric(unlist(tmp))

# focus on intra-tumour sample correlations
for (m in c("t", "b")) {
  distt <- mget(ls(pattern = paste0("dist_m", m)))[[1]]
  tmp <- as.data.frame(as.matrix(distt))
  tmp <- tmp[!grepl("Normal", rownames(tmp)), !grepl("Normal", colnames(tmp))]

  # get column ids for vectors to compare
  seq_i <- seq_len(ncol(tmp))
  seq_j <- seq_len(ncol(tmp))

  # loop over columns pairs ij where j>i
  for (i in seq_i) {
    name_i <- substr(rownames(tmp)[i], 1, 6)
    for (j in seq_j) {
      #extract vector pair
      name_j <- substr(colnames(tmp)[j], 1, 6)
      if (name_i != name_j) {
        # get correlation
        tmp[j, i] <- NA
      }
    }
  }
  vec <- as.numeric(na.omit(unlist(tmp)))
  assign(paste0("dist_", m, m, "_intra"), vec)
}
rm(tmp, vec, distt, seq_i, seq_j, i, j, name_i, name_j)

# focus on inter-tumour sample correlations
for (m in c("t", "b")) {
  distt <- mget(ls(pattern = paste0("dist_m", m)))[[1]]
  tmp <- as.data.frame(as.matrix(distt))
  tmp <- tmp[!grepl("Normal", rownames(tmp)), !grepl("Normal", colnames(tmp))]

  # get column ids for vectors to compare
  seq_i <- seq_len(ncol(tmp))
  seq_j <- seq_len(ncol(tmp))

  # loop over columns pairs ij where j>i
  for (i in seq_i) {
    name_i <- substr(rownames(tmp)[i], 1, 6)
    for (j in seq_j) {
      #extract vector pair
      name_j <- substr(colnames(tmp)[j], 1, 6)
      if (name_i == name_j) {
        # get correlation
        tmp[j, i] <- NA
      }
    }
  }
  vec <- as.numeric(na.omit(unlist(tmp)))
  assign(paste0("dist_", m, m, "_inter"), vec)
}
rm(tmp, vec, distt, seq_i, seq_j, i, j, name_i, name_j)

# focus on tumour-normal
for (m in c("t", "b")) {
  distt <- mget(ls(pattern = paste0("dist_m", m)))[[1]]
  tmp <- as.data.frame(as.matrix(distt))
  tmp <- tmp[grepl("Normal", rownames(tmp)), !grepl("Normal", colnames(tmp))]
  vec <- as.numeric(na.omit(unlist(tmp)))
  assign(paste0("dist_", m, "n"), vec)
}
rm(tmp, vec, distt, m)

# # LUAD vs LUSC
# LUAD <- as.character(patient_info[Histology=="Invasive adenocarcinoma", unique(Patient)])
# LUSC <- as.character(patient_info[Histology!="Invasive adenocarcinoma", unique(Patient)])
# for(m in c("t", "b")){
#   distt <- mget(ls(pattern = paste0("dist_m", m)))[[1]]
#   tmp <- as.data.frame(as.matrix(distt))
#   tmp <- tmp[!grepl("Normal", rownames(tmp)),!grepl("Normal", colnames(tmp))]
#   tmp <- tmp[substr(rownames(tmp), 1,6)%in%LUAD,substr(rownames(tmp), 1,6)%in%LUSC]
#   vec <- as.numeric(na.omit(unlist(tmp)))
#   assign(paste0("dist_", m, m, "_hist"), vec)
# }
# rm(tmp, vec, distt, LUAD, LUSC)

patt <- ls(pattern = "dist")
lst <- mget(patt[!patt %in% c("dist_mb", "dist_mt", "dist_nn_fm", "dist_mtf")])
new_lst <- list()
for (i in seq_len(length(lst))) {
  nam <- names(lst)[i]
  nam <- substr(nam, 6, nchar(nam))
  #extract distances
  x <- as.numeric(lst[[i]])
  #store in nested list of dfs
  flag <- ifelse(grepl("t(t|n)_", nam), "CAMDAC", "bulk")
  new_lst[[i]] <-
    data.frame(x = x,
      flag = rep(nam, length(x)),
      flag2 = rep(flag, length(x)))
}
# bind nested lists items
dt_list <- rbindlist(new_lst)
rm(nam, x, i, patt, lst, new_lst)

# Add the FACS data
tmp <- dist_mtf[, .SD, .SDcols = 2:4]
dt_list <- rbind(dt_list, tmp)

# order flag
flags <- dt_list[, as.character(unique(flag))]
ordr <- c(flags[grepl("^n", flags)],
  "bn", "tn", "bb_intra",   "tt_intra", "bb_inter", "tt_inter")
dt_list[, flag := factor(flag, levels = ordr, ordered = TRUE)]
dt_list <- dt_list[order(flag), ]
rm(flags, ordr)

# annotate sequencing and cancer type
dt_list[, "platform" := "rrbs"]
dt_list[, "cancer_type" := "lung"]

################################################################################
# Load WGBS - mCRPC
################################################################################

# load pairwaise pearson correlation data
dt_wgbs_mcrpc_nn <- fread("inputs/wgbs/mcrpc/normal_pdist.tsv")
dt_wgbs_mcrpc_bn <- fread("inputs/wgbs/mcrpc/bulk_pdist.tsv")
dt_wgbs_mcrpc_tn <- fread("inputs/wgbs/mcrpc/pure_pdist.tsv")
dt_wgbs_mcrpc_bb_inter <- fread("inputs/wgbs/mcrpc/inter_bulk_pdist.tsv")
dt_wgbs_mcrpc_tt_inter <- fread("inputs/wgbs/mcrpc/inter_pure_pdist.tsv")

# combine the data
dt_wgbs_mcrpc <- do.call(rbind, mget(ls(pattern = "dt_wgbs_mcrpc")))

# annotate sequencing and cancer type
dt_wgbs_mcrpc[, "platform" := "wgbs"]
dt_wgbs_mcrpc[, "cancer_type" := "mcrpc"]

################################################################################
# Load WGBS - CRC
################################################################################

# load pairwaise pearson correlation data
dt_wgbs_crc_nn <- fread("inputs/wgbs/crc/normal_pdist.tsv")
dt_wgbs_crc_bn <- fread("inputs/wgbs/crc/bulk_pdist.tsv")
dt_wgbs_crc_tn <- fread("inputs/wgbs/crc/pure_pdist.tsv")
dt_wgbs_crc_bb_inter <- fread("inputs/wgbs/crc/inter_bulk_pdist.tsv")
dt_wgbs_crc_tt_inter <- fread("inputs/wgbs/crc/inter_pure_pdist.tsv")

# combine the data
dt_wgbs_crc <- do.call(rbind, mget(ls(pattern = "dt_wgbs_crc")))

# annotate sequencing and cancer type
dt_wgbs_crc[, "platform" := "wgbs"]
dt_wgbs_crc[, "cancer_type" := "crc"]

################################################################################
# Load RRBS - CRC
################################################################################

# load pairwaise pearson correlation data
dt_rrbs_crc_nn <- fread("inputs/rrbs/crc/normal_pdist.tsv")
dt_rrbs_crc_bn <- fread("inputs/rrbs/crc/bulk_pdist.tsv")
dt_rrbs_crc_tn <- fread("inputs/rrbs/crc/pure_pdist.tsv")
dt_rrbs_crc_bb_inter <- fread("inputs/rrbs/crc/inter_bulk_pdist.tsv")
dt_rrbs_crc_tt_inter <- fread("inputs/rrbs/crc/inter_pure_pdist.tsv")

# combine the data
dt_rrbs_crc <- do.call(rbind, mget(ls(pattern = "dt_rrbs_crc")))

# annotate sequencing and cancer type
dt_rrbs_crc[, "platform" := "rrbs"]
dt_rrbs_crc[, "cancer_type" := "crc"]

################################################################################
# plot the results
################################################################################

# bind all of the data
dt_tmp <- rbind(dt_wgbs_mcrpc, dt_rrbs_crc, dt_wgbs_crc)
dt_new_data <- dt_tmp[, .(x = pdist,
  flag = ifelse(analysis == "bulk_vs_normal", "bn",
         ifelse(analysis == "inter-bulk_vs_bulk", "bb_inter",
         ifelse(analysis == "inter-pure_vs_pure", "tt_inter",
         ifelse(analysis == "normal_vs_normal", "nn",
         ifelse(analysis == "pure_vs_normal", "tn",
         NA))))),
  flag2 = ifelse(grepl("pure", analysis) == TRUE, "CAMDAC", "bulk"),
  platform, cancer_type)]
dt_final <- rbind(dt_list, dt_new_data)
rm(dt_tmp)

# plot ggridges
p <- ggplot(data = dt_final, aes(fill = flag, col = flag, group = flag)) +
  geom_density_ridges(data = dt_final,
                      aes(x = x, y = flag, fill = flag, col = flag,
                          group = flag, height = ..density..),
                      alpha = 0.5, scale = 5) +
  scale_color_manual(name = "sample\npair",
    values = c("nn" = "forestgreen", #"nn_fm"="khaki2",
               "bb_inter" = "grey15", "bb_intra" = "grey55", "bn" = "grey85",
               "tn" = "darkblue","tt_inter" = "dodgerblue", "tt_intra" = "lightblue")) +
  scale_fill_manual(name = "sample\npair",
    values=c("nn" = "forestgreen", #"nn_fm"="khaki2",
             "bb_inter" = "grey15", "bb_intra" = "grey55", "bn" = "grey85",
             "tn" = "darkblue", "tt_inter" = "dodgerblue", "tt_intra" = "lightblue"))+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme_ridges(font_size = 7) +
  theme(legend.text = element_text(size = 6.5),
        axis.text.x = element_text(size = 6.5),
        text = element_text(size = 6.5),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7))

# save grid
f_nm <- "outputs/pearson_corr_distances_bulk_vs_pure_tumour_all_patients_all_samples.pdf"
ggsave(file = f_nm, plot = p, device = "pdf", width = 3.5, height = 3.75, units = "in")

# plot boxplot
p <- ggplot(data = dt_final, aes(fill = flag, col = flag)) +
  geom_boxplot(data = dt_final,
               aes(x = x, y = flag, fill = flag, col = flag),
               alpha = 0.5, outlier.shape = NA) +
  geom_point(data = dt_final[flag2 == "FACS"],
             aes(x = x, y = flag, fill = flag, col = flag, shape = flag2),
             alpha = 0.75, size = 2, position = position_jitterdodge(5)) +
  scale_color_manual(name = "sample\npair",
    values = c("nn" = "forestgreen", #"nn_fm"="khaki2",
               "bb_inter" = "grey15", "bb_intra" = "grey55", "bn" = "grey85",
               "tn" = "darkblue","tt_inter" = "dodgerblue", "tt_intra" = "lightblue")) +
  scale_fill_manual(name = "sample\npair",
    values=c("nn" = "forestgreen", #"nn_fm"="khaki2",
             "bb_inter" = "grey15", "bb_intra" = "grey55", "bn" = "grey85",
             "tn" = "darkblue", "tt_inter" = "dodgerblue", "tt_intra" = "lightblue")) +
  scale_shape_manual(name = "", values = c("CAMDAC" = 16, "FACS" = 17)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  theme(legend.text = element_text(size = 6.5),
        axis.text.x = element_text(size = 6.5),
        text = element_text(size = 6.5),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7))

# save grid
f_nm <- "outputs/pearson_corr_distances_all_samples_new.pdf"
ggsave(file = f_nm, plot = p, device = "pdf",
       width = 3.5, height = 3.75, units = "in")

# 
tmp <- dist_mtf[, .(x = x[flag == "tt_intra"], y = x[flag == "bb_intra"]), by = .(sample)]
tmp[, t.test(x, y, paired = TRUE)]

# END