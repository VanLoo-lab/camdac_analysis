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
x <- c("Biobase", "GenomicRanges", "Rsamtools", "data.table", "stringr",
  "readr", "plyr", "dplyr", "ggforce", "ggplot2", "gridExtra", "betareg",
  "flexmix", "scales", "mixtools")

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
#' @param x Numerical vector of values from which to estimate
#' the mode
#' @return The modes of a multi-modal distribution
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

################################################################################
# Get OC and UC
################################################################################

# set parent directory
tx_dir <- "~/rdeocd/txMeth/"

# set the working directory
work_dir <- "~/data/CAMDAC_paper/figureS5"
setwd(work_dir)

# set the patient id
patient_id = "CRUK0062"

# load the normal patient data
path_patient <- file.path(tx_dir, patient_id)
load(file.path(path_patient, "Methylation/N/gr_normal_m.RData"))

# convert to a data.frame
tmp <- data.frame(gr_normal_m, stringsAsFactors = FALSE)
rm(gr_normal_m)

# investigate normal methylation rates HDI
ggplot(tmp, aes(y = ..density..)) +
  geom_density(aes(x = m_n),  adjust = 0.75, col = "black") +
  geom_density(aes(x = m_n_low), col = "blue") +
  geom_density(aes(x = m_n_high), col = "red") +
  theme_classic()

ggplot(tmp, aes(y=..count..))+
  geom_histogram(aes(x = m_n), col = "black", bins = 100, alpha = 0) +
  geom_histogram(aes(x = m_n_low), col = "blue", bins = 100, alpha = 0) +
  geom_histogram(aes(x = m_n_high), col = "red", bins = 100, alpha = 0) +
  theme_classic()

# get the median 
median(tmp$cov_n)

# get the highest density bins
x <- cut(tmp$m_n, breaks = 50)
head(sort(table(x), decreasing=TRUE))
rm(x)

# get the density
d <- density(tmp$m_n, adjust = 1)

# find the modes
mymodes_indices <- find_modes(d$y)
mymodes <- d$x[mymodes_indices]

# approx shift due to OC and demythylation with age
OC <- mymodes[which(mymodes > 0.85 & mymodes < 0.985)]

# approx shift due to UC
UC <- mymodes[which(mymodes > 0 & mymodes < 0.04)]

suppressWarnings(rm(d, tmp, mymodes_indices, mymodes, find_modes))

################################################################################
# Get bulk methylation plots with CN and purity and DMPs - TRACERx
################################################################################

# set vars
project <- "tracerx"
build <- "hg19"
f_nm <-
  file.path(tx_dir, 
    "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
  )

# load patient_info
patient_info <- fread(f_nm, header = TRUE, stringsAsFactors = FALSE)
rm(f_nm)

# remove lymph nodes
patient_info <- patient_info[!grepl("LN", Sample),]

# get patient ids
patient_ids <- unique(patient_info$Patient)

#
multi_list_patient <- list()
reseg <- TRUE
for(patient_id in patient_ids){
  # set the patient path
  path_patient <- file.path(tx_dir, patient_id, "/")
  
  # make plots
  sample_ids <- unique(patient_info$Sample)[patient_info$Patient == patient_id]
  load(file.path(path_patient,"Methylation/N/gr_normal_m.RData"))
  
  # get normal methylations density
  d <- density(gr_normal_m$m_n, adjust=0.95)
  rm(gr_normal_m)
  
  # find mods
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to OC and demythylation with age
  OC <- mymodes[which(mymodes>0.85&mymodes<0.98)]
  
  # approx shift due to UC
  UC <- mymodes[which(mymodes>0&mymodes<0.04)]
  rm(d, mymodes, mymodes_indices)
  
  multi_list_sample <- list()
  for(sample_id in sample_ids){
    # Load the CN and normal annotated tumour data
    if(reseg==FALSE|patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){ 
      load(
        paste0(path_patient, "/Methylation/", sample_id, "/purified_tumour_DMP_calls.RData")
          )
      }
    if(reseg==TRUE&!patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){ 
      load(
        paste0(path_patient, "/Methylation/", sample_id, "/re_segmented/purified_tumour_DMP_calls.RData")
        )
    }

    # make data.table
    df_mt <- data.table(df_mt)
    
    # get multiplicitity of copy number segments outside X
    df_top_multi <- data.table(sort(df_mt[chromosome != "X", table(multi)], decreasing = TRUE))
    colnames(df_top_multi) <- c("multip", "counts")
    
    # get DMPs at each multiplicitities
    tmp <- df_mt[,.(hyper_counts = sum(!is.na(DMP) & DMP == "hyper" & m_n_high < 0.2), 
                    hypo_counts = sum(!is.na(DMP) & DMP == "hypo" & m_n_low > 0.8)), 
                    by = multi]
    
    # annotate DMP counts
    df_top_multi[, c("hyper_counts", "hypo_counts") := tmp[match(df_top_multi$multip, tmp$multi),2:3]]
    # keep common states (we will be looking at CpGs which are fully 0 or 1 in m_n and looking at DM so need high initial numbers)
    df_top_multi <- df_top_multi[counts>1E5 & hypo_counts>500 & hypo_counts>500,]
    
    #  add results columns
    df_top_multi[, c("b0", "expected_b0", "b1", "expected_b1","t0", "expected_t0", "t1", "expected_t1") := 0] 
    
    top_multi <- sort(df_top_multi$multip)
    cat("Top CpG copy number states :", sep="\n")
    cat(top_multi, sep="\n")
    
    # Load the copy number segments file
    if(reseg==FALSE|patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){copy_number_filename = paste(patient_id,sample_id,"ascat.output.RData",sep=".")
    file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/",copy_number_filename)}
    if(reseg==TRUE&!patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){copy_number_filename = paste(patient_id,sample_id,"ascat.output.resegmented.RData",sep=".")
    file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/re_segmented/",copy_number_filename)}
    load(file_CN)
    p = as.numeric(ascat.output$aberrantcellfraction)
    rm(ascat.output)
    
    for(multi in top_multi[!is.na(top_multi)]){
    #for(multi in c("2+0","2+1","2+2")){
      colss<-ifelse(multi=="2+0","dodgerblue",ifelse(multi=="2+1","orange3","firebrick"))
      
      cat("Generateing histograms - ")
      x <- plot_histogram(df_mt=df_mt, p=p, multip = multi, chr=NULL, OC=OC, UC=UC, reseg=reseg, patient_id=patient_id, sample_id=sample_id)
      cat("histograms generated\n")
      print(x)
      
      df_top_multi[multip==multi, c("b0", "expected_b0", "b1", "expected_b1","t0", "expected_t0", "t1", "expected_t1") :=  as.list(x)] 
    }
    
    j = which(sample_ids==sample_id)
    df_top_multi$patient = patient_id
    df_top_multi$sample = sample_id
    multi_list_sample[[j]] <- df_top_multi
  }
  
  i = which(patient_ids==patient_id)
  multi_list_patient[[i]] <- rbindlist(multi_list_sample)
  rm(multi_list_sample, df_mt, multi)
}
 multi_list <- rbindlist(multi_list_patient)
 
 # load patient_info
 patient_info <- as.data.table(read.delim("~/rdeocd/txMeth/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))
 
 # remove lymph nodes
 patient_info <- patient_info[!grepl("LN", Sample),]
 
 tmp <- data.frame(read_csv("~/rdeocd/txMeth/tx100_RRBS_patients.csv"))
 patient_info$ethnicity <- tmp$Ethnicity[match(patient_info$Patien,tmp$SampleID)]
 patient_info$age <- tmp$Age[match(patient_info$Patien,tmp$SampleID)]
 patient_info$smoking <- tmp$Smoking[match(patient_info$Patien,tmp$SampleID)]
 patient_info$pckyrs <- tmp$Pyrs[match(patient_info$Patien,tmp$SampleID)]
 rm(tmp)
 
 #  add clinical info to multiplicity output
 multi_list$rho <- patient_info$rho_reseg[match(multi_list$patient,patient_info$Patient)]
 multi_list$sex <- patient_info$Sex[match(multi_list$patient,patient_info$Patient)]
 multi_list$hist <- patient_info$Histology[match(multi_list$patient,patient_info$Patient)]
 multi_list$pyrs <- patient_info$pckyrs[match(multi_list$patient,patient_info$Patient)]
 multi_list$smoking <- patient_info$smoking[match(multi_list$patient,patient_info$Patient)]
 multi_list$age <- patient_info$age[match(multi_list$patient,patient_info$Patient)]
 multi_list$ethnicity <- patient_info$ethnicity[match(multi_list$patient,patient_info$Patient)]
 
 #  extract allele specific CN and total CN
 multi_list$nA <- as.numeric(substr(multi_list$multip, 1, 1))
 multi_list$nB <- as.numeric(substr(multi_list$multip, 3, 3))
 multi_list$CN <- multi_list$nA+multi_list$nB
 
 tmp <- multi_list[multi_list$nA<=3&multi_list$nB<=2,]
 p1 <- ggplot(tmp, aes(x=t1, y=expected_t1)) + geom_point(shape=16, size=2, alpha=0.5, col="dodgerblue")+
   xlab("Observed") + ylab("Expected") + theme_minimal()
 p2 <- ggplot(tmp, aes(x=t1-expected_t1, y=..density..)) + geom_density(alpha=0.25, col="dodgerblue") +
   xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 p3 <- ggplot(tmp, aes(x=t0, y=expected_t0)) + geom_point(shape=16, size=2, alpha=0.5, col="navyblue")+
   xlab("Observed") + ylab("Expected") + theme_minimal()
 p4 <- ggplot(tmp, aes(x=t0-expected_t0, y=..density..)) + geom_density(alpha=0.25, col="navyblue") +
   xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 p5 <- ggplot(tmp, aes(x=b1, y=expected_b1)) + geom_point(shape=16, size=2, alpha=0.5, col="grey60")+
   xlab("Observed") + ylab("Expected") + theme_minimal()
 p6 <- ggplot(tmp, aes(x=b1-expected_b1, y=..density..)) + geom_density(alpha=0.25, col="grey60") +
   xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 p7 <- ggplot(tmp, aes(x=b0, y=expected_b0)) + geom_point(shape=16, size=2, alpha=0.5, col="grey30") +
   xlab("Observed") + ylab("Expected") + theme_minimal()
 p8 <- ggplot(tmp, aes(x=b0-expected_b0, y=..density..)) + geom_density(alpha=0.25, col="grey30") +
   xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 
 # make grid 
 gglist1 <- list(p1, p2, p3, p4, p5, p6, p7) ; gglist2 <- list(p8) 
 my_layout <- rbind(cbind(1, 2, 3, 4), cbind(5, 6, 7, 8))
 gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
 
 # save grid
 f_nm <- "~/rdeocd/txMeth_analysis/CAMDAC_bulk_and_pure_DM_peaks_observed_vs_expected.pdf"
 ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height=6, units = "in")
 rm(gglist1, gglist2, my_layout, gd, f_nm, p1, p2)
 
 p1 <- ggplot(tmp, aes(x=factor(multip), y=b1-expected_b1)) +
    geom_boxplot(width=0.8, col="grey60")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p2 <- ggplot(tmp, aes(x=factor(multip), y=b0-expected_b0)) +
    geom_boxplot(width=0.8, col="grey30")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p3 <- ggplot(tmp, aes(x=factor(multip), y=t1-expected_t1)) +
    geom_boxplot(width=0.8, col="dodgerblue")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p4 <- ggplot(tmp, aes(x=factor(multip), y=t0-expected_t0)) +
    geom_boxplot(width=0.8, col="navyblue")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p5 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=b1-expected_b1)) +
    geom_boxplot(width=0.8, col="grey60")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 p6 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=b0-expected_b0)) +
    geom_boxplot(width=0.8, col="grey30")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 p7 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=t1-expected_t1)) +
    geom_boxplot(width=0.8, col="dodgerblue")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 p8 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=t0-expected_t0)) +
    geom_boxplot(width=0.8, col="navyblue")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 
 # make grid 
 gglist1 <- list(p1, p2, p3, p4, p5, p6, p7)
 gglist2 <- list(p8) 
 my_layout <- rbind(cbind(1, 2, 3, 4), cbind(5, 6, 7, 8))
 gd <- arrangeGrob(grobs = c(gglist1, gglist2), layout_matrix = my_layout)
 
 # save grid
 f_nm <- "~/rdeocd/txMeth_analysis/CAMDAC_bulk_and_pure_DM_peaks_observed_vs_expected_rho_ASCN.pdf"
 ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height=6, units = "in")
 rm(gglist1, gglist2, my_layout, gd, f_nm, p1, p2)
 
 tmp <- data.frame(observed_b=c(tmp$b0, tmp$b1), 
                   observed_t=c(tmp$t0, tmp$t1),
                   expected_b=c(tmp$expected_b0, tmp$expected_b1), 
                   expected_t=c(tmp$expected_t0, tmp$expected_t1), 
                   rho = c(tmp$rho, tmp$rho),
                   multip = c(tmp$multip, tmp$multip))
 p1 <- ggplot(tmp, aes(x=factor(multip), y=observed_b-expected_b)) +
    geom_boxplot(width=0.8, col="grey60")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p2 <- ggplot(tmp, aes(x=factor(multip), y=observed_t-expected_t)) +
    geom_boxplot(width=0.8, col="dodgerblue")+
    xlab("ASCN") + ylab("Obs - Exp") + theme_minimal()
 p3 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=observed_b-expected_b)) +
    geom_boxplot(width=0.8, col="grey60")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 p4 <- ggplot(tmp, aes(x=factor(round2(rho, digits=1)), y=observed_t-expected_t)) +
    geom_boxplot(width=0.8, col="dodgerblue")+
    xlab("rho") + ylab("Obs - Exp") + theme_minimal()
 p5 <- ggplot(tmp, aes(x=observed_t, y=expected_t)) +
    geom_point(shape=16, size=2, alpha=0.5, col="dodgerblue")+
    xlab("Observed") + ylab("Expected") + theme_minimal()
 p6 <- ggplot(tmp, aes(x=observed_t-expected_t, y=..density..)) +
    geom_density(alpha=0.25, col="dodgerblue") +
    xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 p7 <- ggplot(tmp, aes(x=observed_b, y=expected_b)) +
    geom_point(shape=16, size=2, alpha=0.5, col="grey60")+
    xlab("Observed") + ylab("Expected") + theme_minimal()
 p8 <- ggplot(tmp, aes(x=observed_b-expected_b, y=..density..)) +
    geom_density(alpha=0.25, col="grey60") +
    xlab("Obs - Exp") + ylab("Density") + theme_minimal()
 
 # make grid 
 gglist1 <- list(p1, p2, p3, p4, p5, p6, p7)
 gglist2 <- list(p8) 
 my_layout <- rbind(cbind(1, 2, 3, 4), cbind(5, 6, 7, 8))
 gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
 
 # save grid
 f_nm <- "~/rdeocd/txMeth_analysis/CAMDAC_bulk_and_pure_DM_peaks_observed_vs_expected_hyper_hypo_combined.pdf"
 ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height=6, units = "in")
 rm(gglist1, gglist2, my_layout, gd, f_nm, p1, p2)
 
 table(multi_list[,.(factor(round2(rho, digits=1)), !is.na(multi_list$t0))])
 table(multi_list[,.(factor(round2(rho, digits=1)), !is.na(multi_list$t1))])
 
 
plot_histogram <- function(df_mt, p, multip, chr=NULL, OC, UC, reseg=TRUE, patient_id, sample_id, plot) {
  b0=NA ; t0=NA; b1=NA ; t1=NA ; multi=multip
  #set path and file prefix
  path_patient <- file.path("/home/ecardieux/rdeocd/txMeth/",patient_id,"/")
  if(reseg == FALSE | patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
    f_prefix <-
      paste0(path_patient, "/Methylation/", sample_id, "/", patient_id, "_", sample_id)
  }
  if(reseg == TRUE & !patient_id %in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
    f_prefix <-
      paste0(path_patient, "/Methylation/", sample_id, "/re_segmented/", patient_id, "_", sample_id)
  }
  
  #subset to relevant multiplicity
  df_mt_sub <- df_mt[!is.na(multi) & as.character(multi) == multip,]
  df_mt_sub <- df_mt_sub[chromosome != "X",]
  
  #subset to relevant chromosome
  if(!is.null(chr)){
    df_mt_sub <- df_mt_sub[chromosome == chr,]
    f_prefix  <- paste0(f_prefix , "_chr", chr)
  }
  
  # transfor multi charactor to factor
  df_mt_sub$multi <- factor(df_mt_sub$multi, levels = multip)

  # get CN
  CN = mean(df_mt_sub$CG_CN)
  nA <- substr(multip, 1, 1)
  nB <- substr(multip, 3, 3)
  
  # Get normal copy number at heterozygous SNPs
  CN_norm <- mean(ifelse(!is.na(df_mt_sub$BAF), 1, 2))
  # CN norm is 1 out of 2 at CG-destroying/-forming heterozygous SNPs
  
  # get expected peaks in m_b
  # for hypermethylated CpGs (mn = 0)
  expected_peak_0 <-
    (OC * p * CN + UC * CN_norm * (1 - p)) / (p * CN + CN_norm * (1 - p))
  expected_peak_0_max <-
    (1 * (p + 0.005) * (CN + 0.5) + 0.2 * CN_norm * (1 - (p + 0.005))) /
    ((p + 0.005) * (CN + 0.5) + CN_norm * (1 - (p + 0.005)))
  expected_peak_0_min <-
    (OC * (p - 0.005) * (CN - 0.5) + 0 * CN_norm * (1 - (p - 0.005))) /
    ((p - 0.005) * (CN - 0.5) + CN_norm * (1 - (p - 0.005)))
   
  # for hypomethylated CpGs (mn = 1)
  expected_peak_1 <-
    (UC * p * CN + OC * CN_norm * (1 - p)) /
    (p * CN + CN_norm * (1 - p))
  expected_peak_1_max <-
    (0 * (p + 0.005) * (CN + 0.5) + 0.8 * CN_norm * (1 - (p + 0.005))) /
    ((p + 0.005) * (CN + 0.5) + CN_norm * (1 - (p + 0.005)))
  expected_peak_1_min = (UC*(p-0.005)*(CN-0.5)+1*CN_norm*(1-(p-0.005)))/((p-0.005)*(CN-0.5)+CN_norm*(1-(p-0.005)))
  #expected_peak_1_min = expected_peak_1 + (expected_peak_1-expected_peak_1_max)
   
  df_mt_0 <- df_mt_sub[df_mt_sub$m_n_high<0.2,]
  # y0 <- unname(table(df_mt_0$m_b<(expected_peak_0+0.04)&df_mt_0$m_b>(expected_peak_0-0.04))["TRUE"]) 
  tmp <- df_mt_0[m_b>expected_peak_0-0.05, .SD, .SDcols="m_b"]
  y0 <- 3*as.numeric(sort(table(cut(tmp$m_b, breaks = length(seq(expected_peak_0,1,0.02)))), decreasing=TRUE)[1])
  
  # get normal methylations density
  d<-density(tmp$m_b[tmp$m_b>0.01], adjust=1)
  
  # find local density maxima
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to OC and demythylation with age
  idx <- which(mymodes<expected_peak_0_max+0.2&mymodes>expected_peak_0_min-0.2)
  b0 = mymodes[idx]
  if(length(b0)>1){
    # # get normal methylations density
    # d<-density(tmp$m_b, adjust=2)
    # 
    # # find local density maxima
    # mymodes_indices <- find_modes(d$y)
    # mymodes <- d$x[mymodes_indices]
    # 
    # # approx shift due to OC and demythylation with age
    # b0 = mymodes[idx]
    # if(length(b0)>1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    b0 <- mymodes[which.max(d$y[mymodes_indices])]  
    # }
  }
  #b0 <- round2(b0, digits=3)
  if(length(b0)<1){b0 <- NA}
  rm(tmp, idx )
  
  df_mt_1 <- df_mt_sub[df_mt_sub$m_n_low>0.8,]
  # y1 <- unname(table(df_mt_1$m_b<(expected_peak_1+0.05)&df_mt_1$m_b>(expected_peak_1-0.05))["TRUE"])
  tmp <- df_mt_1[m_b<expected_peak_1+0.05, .SD, .SDcols="m_b"]
  y1 <- 3*as.numeric(sort(table(cut(tmp$m_b, breaks = length(seq(0,expected_peak_1,0.02)))), decreasing=TRUE)[1])
  
  # get normal methylations density
  d<-density(tmp$m_b[tmp$m_b<0.98], adjust=0.8)
  
  # find mods
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to UC
  idx <- which(mymodes>expected_peak_1_max-0.2&mymodes<expected_peak_1_min+0.2)
  b1 = mymodes[idx]
  if(length(b1)>1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    b1 <- mymodes[which.max(d$y[mymodes_indices])]  
  }
  #b1 <- round2(b1, digits=3)
  if(length(b1)<1){b1 <- NA}
  rm(tmp, d, mymodes, mymodes_indices)

  colss<-ifelse(multip%in%c("2+0","4+0"),"dodgerblue",
                ifelse(multip%in%c("2+1","4+2"),"orange3",
                       ifelse(multip%in%c("1+1","2+2"),"firebrick", 
                              ifelse(multip=="1+0", "grey65", "forestgreen"))))

  if(plot==TRUE){
    p0 <- ggplot(df_mt_0) + geom_histogram(aes(x=m_b, y=..count..),fill=colss,color=NA,binwidth=0.02,alpha=0.5) +
      theme_classic(base_size = 7) + coord_cartesian(ylim=c(0,y0),xlim=c(0,1))+
      #ggtitle(paste(patient_id, sample_id)) +
      #geom_text(data=df_mt_0[1,],aes(x = expected_peak_0+0.04,y=y0,label=round(expected_peak_0,digits=3))) +
      geom_vline(aes(xintercept=expected_peak_0), size = 0.4) +
      geom_vline(aes(xintercept=expected_peak_0_max), linetype="dashed", size = 0.25) +
      geom_vline(aes(xintercept=expected_peak_0_min), linetype="dashed", size = 0.25)
    #df <- ggplot_build(p0)$data[[1]][,c("y","x","PANEL")]
    ggsave(file = paste0(f_prefix, "_bulk_methylation_hyper_multi_", nA, "_", nB,".pdf"),
           plot=p0, device="pdf", width=1.5, height=1, units = "in")
  
    p1 <- ggplot(df_mt_1) + geom_histogram(aes(x=m_b, y=..count..),fill=colss,color=NA,binwidth=0.02,alpha=0.5) +
      theme_classic(base_size = 7) + coord_cartesian(ylim=c(0,y1),xlim=c(0,1))+
      #ggtitle(paste(patient_id, sample_id))+
      geom_vline(aes(xintercept=expected_peak_1), size = 0.4) +
      geom_vline(aes(xintercept=expected_peak_1_max), linetype="dashed", size = 0.25) +
      geom_vline(aes(xintercept=expected_peak_1_min), linetype="dashed", size = 0.25)
    ggsave(file = paste0(f_prefix, "_bulk_methylation_hypo_multi_", nA, "_", nB,".pdf"),
           plot=p1, device="pdf", width=1.5, height=1, units = "in")
  }
  
  # get expected peaks in m_t
  # for hypermethylated CpGs (mn=0)
  peak_0 = (expected_peak_0*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*UC))/(p*CN)
  peak_0_max =  (expected_peak_0_max*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*UC))/(p*CN)
  peak_0_min =  (expected_peak_0_min*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*UC))/(p*CN)
  
  # find local counts maxima
  tmp <- df_mt_0[m_t>0.4, .SD, .SDcols="m_t"]
  y0 <- 3*as.numeric(sort(table(cut(tmp$m_t, breaks = length(seq(0.4,peak_0_max,0.02)))), decreasing=TRUE)[1])
  
  # get normal methylations density
  d<-density(tmp$m_t, adjust=1)
  
  # find local density maxima
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to OC and demythylation with age
  idx <- which(mymodes>expected_peak_0&mymodes>0.65&mymodes<1.25)
  t0 = mymodes[idx]
  if(length(t0)>1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    
    max_y <- which.max(d$y[mymodes_indices])
    t0 <- mymodes[max_y]
  }
  #t0 <- round2(t0, digits=3)
  if(length(t0)<1){t0 <- NA}
  rm(tmp, idx, d, mymodes, mymodes_indices)
  
  # for hypomethylated CpGs (mn=1)
  peak_1 = (expected_peak_1*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*OC))/(p*CN)
  peak_1_max = (expected_peak_1_max*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*OC))/(p*CN)
  peak_1_min = (expected_peak_1_min*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*OC))/(p*CN)
  
  tmp <- df_mt_1[m_t<0.6, .SD, .SDcols="m_t"]
  y1 <- 3*as.numeric(sort(table(cut(tmp$m_t, breaks = length(seq(peak_1_max,0.7,0.02)))), decreasing=TRUE)[1])
  
  # get normal methylations density
  d<-density(tmp$m_t, adjust=1)
  
  # find mods
  mymodes_indices <- find_modes(d$y)
  mymodes <- d$x[mymodes_indices]
  
  # approx shift due to UC
  idx <- which(mymodes>-0.25&mymodes<0.25)
  t1 = mymodes[idx]
  if(length(t1)>1){
    mymodes_indices <- mymodes_indices[idx]
    mymodes <- mymodes[idx]
    t1 <- mymodes[which.max(d$y[mymodes_indices])]
  }
  #t1 <- round2(t1, digits=3)
  if(length(t1)<1){t1 <- NA}
  rm(tmp, idx, d, mymodes, mymodes_indices)
  
  if(plot==TRUE){
  p0 <- ggplot(df_mt_0) + geom_histogram(aes(x=m_t, y=..count..),fill=colss,color=NA,binwidth=0.02,alpha=0.5) +
    theme_classic(base_size = 7) + coord_cartesian(ylim=c(0,y0),xlim=c(0,peak_0_max+0.05))+
    #ggtitle(paste(patient_id, sample_id)) +
    geom_vline(aes(xintercept=peak_0), size = 0.4) +
    geom_vline(aes(xintercept=peak_0_max), linetype="dashed", size = 0.25) +
    geom_vline(aes(xintercept=peak_0_min), linetype="dashed", size = 0.25)
    ggsave(file = paste0(f_prefix, "_pure_methylation_hyper_multi_", nA, "_", nB,".pdf", sep = ""),
           plot=p0, device="pdf", width=1.5, height=1, units = "in")

  p1 <- ggplot(df_mt_1) + geom_histogram(aes(x=m_t, y=..count..),fill=colss,color=NA,binwidth=0.02,alpha=0.5) +
    theme_classic(base_size = 7) + coord_cartesian(ylim=c(0,y1),xlim=c(peak_1_max-0.05,1))+
    #ggtitle(paste(patient_id, sample_id)) +
    geom_vline(aes(xintercept=peak_1), size = 0.4) +
    geom_vline(aes(xintercept=peak_1_max), linetype="dashed", size = 0.25) +
    geom_vline(aes(xintercept=peak_1_min), linetype="dashed", size = 0.25)
    ggsave(file = paste0(f_prefix, "_pure_methylation_hypo_multi_", nA, "_", nB,".pdf", sep = ""),
           plot=p1, device="pdf", width=1.5, height=1, units = "in")
  }
  
  return(c(b0=b0, expected_b0=expected_peak_0, b1=b1, expected_b1=expected_peak_1,
           t0=t0, expected_t0=peak_0, t1=t1, expected_t1=peak_1))
}

#p0 + facet_zoom(ylim = c(0, y0))

# df_mt_1 <- df_mt_1[!is.na(df_mt_1$prob_hypo_corr),]
# df_mt_1 <- df_mt_1[df_mt_1$prob_hypo_corr<0.01,]
# p1 <- ggplot(df_mt_1) + geom_histogram(aes(x=m_b, y=..count..,fill=multi,color=multi),bins=50,alpha=0.25) +
#   facet_grid(multi ~ .,scales = "free") + theme_classic()+ggtitle(paste(patient_id, sample_id))+
#   geom_vline(aes(xintercept=expected_peak_1), linetype="dashed")+scale_x_reverse()+coord_cartesian(xlim=c(0,1))
# ggsave(file = paste(path_patient,"Allelecounts/",sample_id,"/Methylation/",patient_id, "_", sample_id, "_bulk_Methylation_multi_tumour_hypo_DMPs_", multi,".pdf", sep = ""),plot=p1,device="pdf", width=3,height=1.5, units = "in")
# 
# df_mt_0 <- df_mt_0[!is.na(df_mt_0$prob_hyper_corr),]
# df_mt_0 <- df_mt_0[df_mt_0$prob_hyper_corr<0.05,]
# p0 <- ggplot(df_mt_0) + geom_histogram(aes(x=m_b, y=..count..,fill=multi,color=multi),bins=50,alpha=0.25) +
#   facet_grid(multi ~ .,scales = "free") + theme_classic()+ggtitle(paste(patient_id, sample_id)) + coord_cartesian(xlim=c(0,1))+
#   geom_vline(aes(xintercept=expected_peak_0), linetype="dashed")
# ggsave(file = paste(path_patient,"Allelecounts/",sample_id,"/Methylation/",patient_id, "_", sample_id, "_bulk_Methylation_multi_tumour_hyper_DMPs_", multi,".pdf", sep = ""),plot=p0,device="pdf", width=3,height=1.5, units = "in")


patient_files <- list.dirs("~/rdeocd/txMeth/",recursive = FALSE)
patient_ids <- unlist(lapply(1:length(patient_files ), function(i, patient_files ) strsplit(patient_files , split="[//]")[[i]][7], patient_files ));rm(patient_files )
patient_ids <- patient_ids[patient_ids!=resolve_tx_patient_id("CRUK0047", .tx_map)]
 
for(patient_id in patient_ids){
path_patient= file.path("/home/ecardieux/rdeocd/txMeth/",patient_id,"/")
load(file.path(path_patient,"Methylation/N/gr_normal_m.RData"))
d<-density(gr_normal_m$m_n, adjust=0.9)
rm(gr_normal_m)
if(patient_id == resolve_tx_patient_id("CRUK0034", .tx_map)){ 
  png("~/rdeocd/txMeth_analysis/all_normal_densities.png")
  plot(d, xlim=c(0,1), ylim=c(0,10))
      } else {
  lines(d)
      }
# find mods
mymodes_indices <- find_modes(d$y)
mymodes <- d$x[mymodes_indices]

# approx shift due to OC and demythylation with age
OC = mymodes[which(mymodes>0.85&mymodes<0.985)]
abline(v=OC, lty="dashed", lwd=0.5)

if(patient_id == resolve_tx_patient_id("CRUK0014", .tx_map)){
    dev.off()
}
}

###############################################
##            Compare M/UM tails             ##
###############################################

# add additional info
tmp <- data.frame(read_csv("~/rdeocd/txMeth/tx100_RRBS_patients.csv"))
tmp <- tmp[tmp$SampleID!="CRUK0047",]

patients1 <- unique(tmp$SampleID)
list_patients_d <- list()
for(patient1 in patients1){
  i <- which(patients1 == patient1)
  set.seed(1234)
    
  # set patient_id2
  tmp2 <- tmp[tmp$SampleID==patient1,]
  
  # set patient_id2
  patients2 <- tmp[tmp$SampleID!=patient1 &
                  tmp$Gender==tmp2$Gender &
                  tmp$Smoking==tmp2$Smoking,"SampleID"]
  patients2 <- patients2[patients2!="CRUK0047"]
  
  list_samples_d <- list()
  for(patient2 in patients2){
    #patient2 <- sample(patient2[patient_id != resolve_tx_patient_id("CRUK0047", .tx_map)], 1, replace=TRUE)
    
    # set patient paths
    path_patient1<-paste0("/home/ecardieux/rdeocd/txMeth/", patient1, "/",sep="")
    path_patient2<-paste0("/home/ecardieux/rdeocd/txMeth/", patient2, "/",sep="")
    
    # Get purity sample 1
    if(reseg==FALSE|patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
       file_n1 = paste0(path_patient1, "Methylation/N/gr_normal_m.RData")
       file_n2 = paste0(path_patient2, "Methylation/N/gr_normal_m.RData")
       # file_DMP1 = paste0(path_patient1, "Methylation/",sample_id, 
       #                   "/purified_tumour_DMP_calls.RData")
       # file_DMP2 = paste0(path_patient2, "Methylation/",sample_id, 
       #                   "/purified_tumour_DMP_calls.RData")
    }
    if(reseg==TRUE&!patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
       file_n1 = paste0(path_patient1, "Methylation/N/gr_normal_m.RData")
       file_n2 = paste0(path_patient2, "Methylation/N/gr_normal_m.RData")
       # file_DMP1 = paste0(path_patient1, "Methylation/",sample_id, 
       #                   "/re_segmented/purified_tumour_DMP_calls.RData")
       # file_DMP2 = paste0(path_patient2, "Methylation/",sample_id, 
       #                   "/re_segmented/purified_tumour_DMP_calls.RData")
    }
    
    # load tumour sample 1
    load(file_n1)
    gr_normal_m1 <- gr_normal_m[as.character(seqnames(gr_normal_m))!="X"]
    
    # load tumour sample 2
    load(file_n2)
    gr_normal_m2 <- gr_normal_m[as.character(seqnames(gr_normal_m))!="X"]
    rm(file_n1, file_n2, gr_normal_m)
    
    gr_normal_m <- mergeByOverlaps(gr_normal_m1, gr_normal_m2)
    gr_normal_m <- gr_normal_m[gr_normal_m$gr_normal_m1$m_n_high<0.2|
                               gr_normal_m$gr_normal_m1$m_n_low>0.8,]
    
    ## data.frame
    d1 <- data.frame(x = c(gr_normal_m$gr_normal_m1$m_n, 
                          gr_normal_m$gr_normal_m2$m_n), 
                    flag = c(rep("matched", length(gr_normal_m$gr_normal_m1)), 
                             rep("unmatched", length(gr_normal_m$gr_normal_m2))), 
                    stringsAsFactors = FALSE)
    d2 <- d1
    rm(gr_normal_m1, gr_normal_m2)
    
    ## separate histograms for both matched and unmatched data
    pb <-  ggplot(d1) + theme_minimal() + theme(legend.position = "none") + 
      geom_histogram(aes(x=x, y=..count.., fill=flag), 
                     binwidth = 0.01, alpha=0.5, position="identity") + 
      scale_fill_manual(name="", values=c("matched"=hcl(0, 120, 80),
                                          "unmatched"=hcl(240, 120, 80))) 
    
    Sys.sleep(1)
    pb_zoom <- ggplot(d2) + theme_minimal() + 
      geom_histogram(aes(x=x, y=..count.., fill=flag), 
                     binwidth = 0.01, alpha=0.5, position="identity") + 
      scale_fill_manual(name="", values=c("matched"=hcl(0, 120, 80),
                                          "unmatched"=hcl(240, 120, 80))) +
      coord_cartesian(ylim=c(0,3E3))
    
    Sys.sleep(1)
    
    gglist1 <- list(pb)
    gglist2 <- list(pb_zoom)
    my_layout <- cbind(1,1,1,2,2,2,2)
    gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
    
    f_nm <- paste0(path_patient1, "Methylation/N/UM_and_M_tail_matched_vs_unmatched_normal_", patient2 ,".pdf")
    ggsave(file = f_nm, plot=gd, device = "pdf", width=8, height = 3.5, units = "in")
    
    # extract m_n unmatched for modelling
    x0 <- gr_normal_m[gr_normal_m$gr_normal_m1$m_n_high<0.2,]
    x <- x0$gr_normal_m1$m_n
    x <- na.omit(x)
    
    # # model beta
    # m1 <- betareg(x ~ 1, data = data.frame(x=[x>0&x<1]))
    # mu <- unname(plogis(coef(m1)[1]))
    # phi <- unname(exp(coef(m1)[2])) # bound by [0,Inf]
    
    # d3 <- density(x)
    # d3 <- data.table(x=d3$x, y=d3$y, ID=round2(d3$x, digits = 2))
    # d3 <- d3[, .(x=mean(x), y=mean(y)), by=.(ID)]
    d3 <- data.table(x=x, ID=as.character(round2(x, digits = 2)))
    d3 <- d3[, .(x=mean(x), counts=length(x)), by=.(ID)]
    d3 <- d3[order(as.numeric(ID)),]
    d3 <- d3[d3$x>=0&d3$x<=1]
    # d3$y <- d3$y/sum(d3$y)
    
    # extract m_n unmatched for modelling
    x <- x0$gr_normal_m2$m_n
    x <- na.omit(x)
    
    d4 <- data.table(x=x, ID=as.character(round2(x, digits = 2)))
    d4 <- d4[, .(x=mean(x), counts=length(x)), by=.(ID)]
    d4 <- d4[d4$x>=0&d4$x<=1]
    d4 <- d4[order(as.numeric(ID)),]
    
    # extract m_n unmatched for modelling
    x1 <- gr_normal_m[gr_normal_m$gr_normal_m1$m_n_low>0.8,]
    x <- x1$gr_normal_m1$m_n
    x <- na.omit(x)
    
    d5 <- data.table(x=x, ID=as.character(round2(x, digits = 2)))
    d5 <- d5[, .(x=mean(x), counts=length(x)), by=.(ID)]
    d5 <- d5[d5$x>=0&d5$x<=1]
    d5 <- d5[order(as.numeric(ID)),]
    
    # extract m_n unmatched for modelling
    x <- x1$gr_normal_m2$m_n
    x <- na.omit(x)
    
    d6 <- data.table(x=x, ID=as.character(round2(x, digits = 2)))
    d6 <- d6[, .(x=mean(x), counts=length(x)), by=.(ID)]
    d6 <- d6[d6$x>=0&d6$x<=1]
    d6 <- d6[order(as.numeric(ID)),]
    rm(x0, x, x1)
    
    df <- data.frame(ID=as.character(seq(0,1,0.01)))
    df$counts_matched_UM <- 0
    df$counts_unmatched_UM <- 0
    df$counts_matched_M <- 0
    df$counts_unmatched_M <- 0
    
    df$counts_matched_UM[!is.na(match(df$ID,d3$ID))] <- d3$counts
    df$counts_unmatched_UM[!is.na(match(df$ID, d4$ID))] <- d4$counts
    df$tail_abs_UM <- df$counts_unmatched_UM - df$counts_matched_UM
    df$tail_perc_UM <- (df$counts_unmatched_UM - df$counts_matched_UM + 1) / (df$counts_matched_UM + 1)
  
    df$counts_matched_M[!is.na(match(df$ID, d5$ID))] <- d5$counts
    df$counts_unmatched_M[!is.na(match(df$ID, d6$ID))] <- d6$counts
    df$tail_abs_M <- df$counts_unmatched_M - df$counts_matched_M
    df$tail_perc_M <- (df$counts_unmatched_M - df$counts_matched_M + 1) / (df$counts_matched_M + 1)
    
    j <- which(patients2==patient2)
    list_samples_d[[j]] <- df
  }
  list_patients_d[[i]] <- rbindlist(list_samples_d)
}

df_summary <- rbindlist(list_patients_d)
df_summary[ , ID := as.character(ID)]
df_summary <- df_summary[, lapply(.SD, median), .SDcols=colnames(df_summary)[grepl("counts_", colnames(df_summary))], by=.(ID)]

df_summary[, tail_perc_UM := (counts_unmatched_UM - counts_matched_UM) / (counts_matched_UM + 1E-5)]
df_summary[, tail_perc_M := (counts_unmatched_M - counts_matched_M) / (counts_matched_M + 1E-5)]

df_summary[, tail_abs_UM := (counts_unmatched_UM - counts_matched_UM) ]
df_summary[, tail_abs_M := (counts_unmatched_M - counts_matched_M) ]

df_summary[, d_UM := (counts_unmatched_UM)/sum(counts_unmatched_UM)]
df_summary[, d_M := (counts_unmatched_M)/sum(counts_unmatched_M)]

# df_summary[, tail_perc_UM := ifelse(tail_perc_UM<0, 0, tail_perc_UM)]
# df_summary[, tail_perc_M := ifelse(tail_perc_M<0, 0, tail_perc_M)]

## separate histograms for both matched and unmatched data
p1 <-  ggplot(df_summary) + 
  theme_minimal() + theme(legend.position = "none") + ylab("log2 fold change") +
  geom_col(aes(x=as.numeric(ID), y=log2(abs(tail_perc_UM))), fill=hcl(240, 120, 80), alpha=0.5)

p2 <-  ggplot(df_summary) + 
  theme_minimal() + theme(legend.position = "none") + ylab("log2 fold change") +
  geom_col(aes(x=as.numeric(ID), y=log2(abs(tail_perc_M))), fill=hcl(0, 120, 80), alpha=0.5) 

p3 <-  ggplot(df_summary) + 
  theme_minimal() + theme(legend.position = "none") + 
  geom_col(aes(x=as.numeric(ID), y=counts_unmatched_UM), fill=hcl(240, 120, 80), alpha=0.5) +
  geom_col(aes(x=as.numeric(ID), y=counts_unmatched_M), fill=hcl(0, 120, 80), alpha=0.5) 

Sys.sleep(1)
p3_zoom <- ggplot(df_summary) + 
  theme_minimal() + theme(legend.position = "none") + 
  geom_col(aes(x=as.numeric(ID), y=counts_unmatched_UM), fill=hcl(240, 120, 80), alpha=0.5) +
  geom_col(aes(x=as.numeric(ID), y=counts_unmatched_M), fill=hcl(0, 120, 80), alpha=0.5) +
  coord_cartesian(ylim=c(0,3E3))

Sys.sleep(1)

gglist1 <- list(p1, p2, p3)
gglist2 <- list(p3_zoom)
my_layout <- rbind(cbind(1,2),cbind(3,4))
gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)

f_nm <- paste0("~/rdeocd/txMeth_analysis/UM_and_M_tail_matched_vs_unmatched_normal.pdf")
ggsave(file = f_nm, plot=gd, device = "pdf", width=8, height = 5, units = "in")

rm(tmp, tmp2, d, d1, d2, d3, d4, d5, d6, gr_normal_m, pb, pb_zoom, 
   gglist1, gglist2, gd, my_layout, f_nm, p1, p2,  p3_zoom)
rm(list=ls(pattern="patients"))

###############################################
##             SIMULATION v2                 ##
###############################################

# Calculate HDI by simulation
#
# Computes highest density interval from a sample of representative values,
#   estimated as shortest credible interval for a unimodal distribution
# Arguments:
#' @param M_b vector with counts methylated in the tumour
#' @param UM_b vector with counts unmethylated in the tumour
#' @param M_n vector with counts methylated in the normal
#' @param UM_n vector with counts unmethylated in the normal
#' @param p tumour purity
#' @param CN vector with total copy number for each tumour CpG loci
#' @param credMass default is 0.99
#' is a scalar between 0 and 1, indicating the mass within the credible interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
HDIofMCMC = function(M_b,UM_b,M_n,UM_n,p,CN,credMass=0.99) {
  
  # sampleVec is a vector of representative values for the methylation rate 
  # probability distribution at a CpG
  sampleVec1 <- rbeta(n = 1000, shape1 = M_b+1, shape2 = UM_b+1) 
  sampleVec2 <- rbeta(n = 1000, shape1 = M_n+1, shape2 = UM_n+1)
  
  # get the mb and mn constant
  cons1 <- p*CN + 2*(1-p)
  cons2 <- 2*(1-p)
  
  # obtain the mt sample vector
  sampleVec_scaled <-  (cons1*sampleVec1 - cons2*sampleVec2)/(p*CN)
  
  # sort the mt sample vector
  sortedPts <- sort(sampleVec_scaled)
  
  # get the width of the nth percentile where n=credMass
  ciIdxInc <- ceiling(credMass * length(sortedPts)) 
  
  # get the diff between all pairs at a suitable width
  ciWidth <- diff(x = sortedPts, lag = ciIdxInc) 
  
  # calculate mt HDI
  HDIlim <- c(lower = sortedPts[which.min(ciWidth)], 
              upper = sortedPts[which.min(ciWidth)+ciIdxInc])
  
  return(HDIlim)
}


# Arguments:
#'  @param ICDFname is R's name for the inverse cumulative density function
#' of the distribution.
#' @param credMass is the desired mass of the HDI region.
#' @param tol is passed to R's optimize function, 
#' the lower the tolerance,the longer the optimisation, but the higher the accuracy.
#' tol=1e-4 gives values of the same accurary as our max resolution
#' Return value:
#' Highest density iterval (HDI) limits in a vector.
#' Example of use: For determining HDI of a beta(30,12) distribution, type
#' HDIofICDF( qbeta , shape1 = 30+1 , shape2 = 12+1 )
#' Notice that the parameters of the ICDFname must be explicitly named;
#' e.g., HDIofICDF( qbeta , 30+1 , 12+1 ) does not work.
#' Adapted and corrected from Greg Snow's TeachingDemos package.

# Source fct outside of loop to speed up code
intervalWidth =  function(lowTailPr,ICDFname,credMass, ... ) {
  ICDFname(credMass+lowTailPr, ... ) - ICDFname(lowTailPr, ... )
}

HDIofICDF = function(ICDFname, credMass=0.99 , tol=1e-4, ... ) {
  
  incredMass = 1.0 - credMass
  
  optInfo = optimize(f = intervalWidth, interval = c(0,incredMass) , ICDFname=ICDFname , credMass=credMass , tol=tol , ... )
  
  HDIlowTailPr = optInfo$minimum
  vec <- setNames(object = ICDFname(c(HDIlowTailPr, credMass+HDIlowTailPr), ... ), nm = c("low", "high"))
  return(vec)
}

overlaps <- function(v, w, x, y) {
  value = FALSE
  if(sum(!is.na(c(v, w, x, y)))==4){
  if(x<=v&v<=y){ value = TRUE }
  if(x<=w&w<=y){ value = TRUE }
  }
  return(value)
}

# get overlaps
# Evan Miller's closed form solution for the probability that a draw from a beta dist is greater than anoter, 
# p_A ~ Beta(alpha_A, beta_A), p_B ~ Beta(alpha_B, beta_B). In this case, A is the normal and B is the bulk tumour. 
# Alpha is the counts methylated, beta is counts unmethylated
h <- function(alpha_n, beta_n, alpha_b, beta_b){
  j <- seq.int(0, round(alpha_b)-1)
  log_vals <- (lbeta(alpha_n + j, beta_n + beta_b) - log(beta_b + j) -
                 lbeta(1 + j, beta_b) - lbeta(alpha_n, beta_n))
  1 - sum(exp(log_vals))
}

# Set static vars
n_cores = 16
set.seed(1234)
# results_list <- list()
counter = 0
reseg=TRUE
dt_top_multi_list = list()

# load all normal CpGs
load("~/storage/CI_FeberTanic_TracerX/results_211120/all_normal_CpGs_matrix.RData")

# get patient_ids
patient_files <- list.dirs("~/rdeocd/txMeth/",recursive = FALSE)
patient_ids <- unlist(lapply(1:length(patient_files ), function(i, patient_files ) strsplit(patient_files , split="[//]")[[i]][7], patient_files ));rm(patient_files )
patient_ids <- patient_ids[patient_ids!=resolve_tx_patient_id("CRUK0047", .tx_map)]

# get patient info
for(patient_id in patient_ids){ 
  path_patient<-paste0("/home/ecardieux/rdeocd/txMeth/", patient_id, "/",sep="")
  setwd(path_patient)
  
  files <- list.dirs(paste0("~/rdeocd/txMeth/",patient_id,"/archive/Allelecounts"),recursive = FALSE)
  sample_ids <- unlist(lapply(1:length(files), function(i, files) strsplit(files, split="[//]")[[i]][9], files));rm(files)
  sample_ids <- sample_ids[grepl(pattern="^R[0-9]", sample_ids)|grepl(pattern="^LN[0-9]", sample_ids)]
  
  n <- length(sample_ids) 
  
  load(file.path(path_patient,"Methylation/N/gr_normal_m.RData"))
  
  # Extract methylated peak from bimodal data
  set.seed(1234)
  d <- data.frame(y = c(sample(gr_normal_m$m_n[gr_normal_m$m_n_low > 0.8 & gr_normal_m$m_n != 1], 
                               1E3, replace=FALSE))) 
 
  # fit beta mixture model
  m <- betamix(y ~ 1 | 1, data = d, k = 1)
  
  # extract mu and phi parameters
  OC <- unname(plogis(coef(m)[1])) # bound by [0,1]

  # Extract unmethylated peak from bimodal data
  d <- data.frame(y = sample(gr_normal_m$m_n[gr_normal_m$m_n_high < 0.2 & gr_normal_m$m_n != 0], 
                               1E3, replace=FALSE)) 
  
  # fit beta mixture model
  m <- betamix(y ~ 1 | 1, data = d, k = 1)
  
  # extract mu and phi parameters
  UC <- unname(plogis(coef(m)[1])) # bound by [0,1]

  # get sample info
  for(sample_id in sample_ids){
    
    # Get purity sampple 1
    if(reseg==FALSE|patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
      copy_number_filename = paste(patient_id, sample_id,"ascat.output.RData",sep=".")
      file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/",copy_number_filename)
      file_DMP = paste0(path_patient, "Methylation/",sample_id, "/purified_tumour_DMP_calls.RData")
    }
    if(reseg==TRUE&!patient_id%in% c(resolve_tx_patient_id("CRUK0010", .tx_map), resolve_tx_patient_id("CRUK0047", .tx_map))){
      copy_number_filename = paste(patient_id, sample_id,"ascat.output.resegmented.RData",sep=".")
      file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/re_segmented/",copy_number_filename)
      file_DMP = paste0(path_patient, "Methylation/",sample_id, "/re_segmented/purified_tumour_DMP_calls.RData")
    }
    load(file_CN) 
    p = as.numeric(ascat.output$aberrantcellfraction)
    rm(ascat.output, file_CN, copy_number_filename)
        
    # load tumour sample
    load(file_DMP)
    
    # select relevant cols
    cols <-  c("chromosome", "start", "end", "cov_b", "CN", "nA", "nB")
    df_mt <- df_mt[,c(cols,"cov_n","m_n","m_n_low","m_n_high", "m_t", "m_t_low",
                      "m_t_high", "m_b", "m_b_low","m_b_high", "multi")]
    df_mt <- df_mt[df_mt$chromosome!="X",]
    rm(file_DMP, cols)
    
    # store copy number
    CN <- df_mt$CN
    nA <- df_mt$nA
    nB <- df_mt$nB
    
    # load all normal CpGs
    flag <- as.logical(all_normal_samples[, !is.na(.SD) , .SDcols=paste0(patient_id, "_N")])
    
    # remove NAs in matched normal
    tmp <- all_normal_samples[flag==TRUE,]  
    
    # remove matched normal and keep all unmatched
    # normals should be more or less the same 
    # this will be used to model normal noise in bulk tumour
    cols <- colnames(tmp)[!grepl(patient_id, colnames(tmp))]
    tmp <- tmp[, .SD, .SDcols = cols]
    
    # get unmatched normal CpGs which are UM in matched normal
    UM_coord <- df_mt[df_mt$m_n_high<0.2, 1:3]
    UM <- merge(tmp, UM_coord)
    
    # get unmatched normal CpGs which are M in matched normal
    M_coord <- df_mt[df_mt$m_n_low>0.8, 1:3]
    M <- merge(tmp, M_coord)
   
    # flag columns with normal methylation information
    cols <- colnames(tmp)[grepl("^[A-Za-z]+[0-9]", colnames(tmp))]
    
    # get unmatched normal
    UM$sampling <- sample(cols, nrow(UM), replace = TRUE)
    UM <- UM[, cbind(.SD, chosen = .SD[[.BY$sampling]]), by=sampling]$chosen
    UM <- sample(UM[!is.na(UM)], length(UM), replace = TRUE)
    # UM <- as.numeric(UM[, rowMeans(.SD, na.rm = TRUE), .SDcols=cols])
    
    M$sampling <- sample(cols, nrow(M), replace = TRUE)
    M <- M[, cbind(.SD, chosen =.SD[[.BY$sampling]]), by=sampling]$chosen
    M <- sample(M[!is.na(M)], length(M), replace = TRUE)
    rm(flag, UM_coord, M_coord, tmp, cols)
    
    # get normal methylation rate vectors for confidently unmethylated and for methylated loci
    # UM <- df_mt$m_n[df_mt$m_n_high<0.2]
    # M  <- df_mt$m_n[df_mt$m_n_low>0.8]

    # Assign a normal methylation state, UM or M, to all CpGs 
    # The selector randomises UM and M ground state respecting the UM to M ratio
    selector <- c(rep(0,length(M)),rep(1,length(UM)))
    selector <- selector[sample(1:length(selector),nrow(df_mt),replace = T)]
        
    # extract normal methylation state from UM and M 
    prop_n <- numeric(length=nrow(df_mt))
    prop_n[selector == 0] <- M[sample(1:length(M),sum(selector==0),replace = T)]
    prop_n[selector == 1] <- UM[sample(1:length(UM),sum(selector==1),replace = T)]
        
    # Assign the tumour cells methylation state, UM or M, to all CpGs
    # This one will be used to simulate CpGs on the major allele (nA)
    selector1 <- c(rep(0,length(M)),rep(1,length(UM)))
    selector1 <- selector1[sample(1:length(selector1),nrow(df_mt),replace = T)]

    # extract the tumour methylation state from UM and M for the major allele
    prop_t1 <- numeric(length=nrow(df_mt))
    prop_t1[selector1 == 0] <- M[sample(1:length(M),sum(selector1==0),replace = T)]
    prop_t1[selector1 == 1] <- UM[sample(1:length(UM),sum(selector1==1),replace = T)]
    
    # Assign the tumour cells methylation state, UM or M, to all CpGs
    # This one will be used to simulate CpGs on nB 
    selector2 <- c(rep(0,length(M)),rep(1,length(UM)))
    selector2 <- selector2[sample(1:length(selector2),nrow(df_mt),replace = T)]
    
    # extract the tumour methylation state from UM and M for the minor allele
    prop_t2 <- numeric(length=nrow(df_mt))
    prop_t2[selector2 == 0] <- M[sample(1:length(M),sum(selector2==0),replace = T)]
    prop_t2[selector2 == 1] <- UM[sample(1:length(UM),sum(selector2==1),replace = T)]
    
    # # Assign the tumour cells methylation state, UM or M, to all CpGs
    # # This one will be used to simulate copy specific changes only
    # selector3 <- c(rep(0,length(M)),rep(1,length(UM)))
    # selector3 <- selector3[sample(1:length(selector3),nrow(df_mt),replace = T)]
    # 
    # # extract the tumour methylation state from UM and M
    # # for the copy-specific epimutations
    # flag <- selector1 == selector2 & selector == selector1 & nA == 2 & nB == 0
    # prop_t3 <- rep(0, nrow(df_mt))
    # prop_t3[selector3 == 0 & flag == TRUE] <- 
    #   M[sample(1:length(M), sum(selector3 == 0 & flag == TRUE),replace = T)]
    # prop_t3[selector3 == 1 & flag == TRUE] <- 
    #   UM[sample(1:length(UM), sum(selector3 == 1 & flag == TRUE),replace = T)]

    # in 2+0, allow for copy-specific methylation
    flag_n <- nA == 2 & nB == 0
    n <- ifelse(flag_n == TRUE, sample(1:2, 1), nA)
    # get cov_t for the major allele
    cov_t1 <- unlist( mcmapply( 
      function (p, CN, cov, n) { 
        prop <- p*n / (p*CN + 2*(1-p)) # if CN=2, prop = purity
        cov_t <- rbinom(1,cov,prop)
        return(cov_t)
      },p=p, CN=df_mt$CN, cov=df_mt$cov_b, n=n, 
      mc.cores = n_cores))  
    
    # get counts methylated in the tumour (methylation prior = normal state)
    M_t1 <- unlist(mcmapply( function (cov, prop) rbinom(1,cov,prop) , 
                             cov=cov_t1, prop=prop_t1, mc.cores = n_cores)) 
    UM_t1 <- cov_t1 - M_t1
    m_t1 <- M_t1 /(M_t1 + UM_t1)
    
    # in 2+0, allow for copy-specific methylation
    n <- ifelse(flag_n == TRUE, sample(0:1, 1), nB)

    # get cov_t for the minor allele
    cov_t2 <- unlist( mcmapply( 
      function (p, CN, cov, n) { 
        prop <- p*n / (p*CN + 2*(1-p)) # if CN=2, prop = purity
        cov_t <- rbinom(1,cov,prop)
        return(cov_t)
      },p=p, CN=df_mt$CN, cov=df_mt$cov_b, n=n, 
      mc.cores = n_cores))  
    
    # get counts methylated in the tumour (methylation prior = normal state)
    M_t2 <- unlist(mcmapply( function (cov, prop) rbinom(1,cov,prop) , 
                             cov=cov_t2, prop=prop_t2, mc.cores = n_cores)) 
    UM_t2 <- cov_t2 - M_t2
  
    # # get cov_t for copy specific changes
    # cov_t3 <- rep(0, nrow(df_mt))
    # cov_t3[flag==TRUE] <- unlist( mcmapply(
    #   function (p, CN, cov, n) {
    #     prop <- p*n / (p*CN + 2*(1-p)) # if CN=2, prop = purity
    #     cov_t <- rbinom(1,cov,prop)
    #     return(cov_t)
    #   },p=p, CN=df_mt$CN[flag==TRUE], cov=df_mt$cov_b[flag==TRUE], n=1,
    # mc.cores = n_cores))
    # 
    # # get counts methylated in the tumour (methylation prior = normal state)
    # M_t3 <- rep(0, nrow(df_mt))
    # M_t3[flag==TRUE] <- unlist(mcmapply( function (cov, prop) rbinom(1,cov,prop) ,
    #                          cov=cov_t3[flag==TRUE], prop=prop_t3[flag==TRUE], mc.cores = n_cores))
    # UM_t3 <- cov_t3 - M_t3

    # get total counts methylated in the tumour
    M_t <- M_t1 + M_t2
    #M_t[flag==TRUE] <- round2(M_t1[flag==TRUE]/2, digits=0) + M_t3[flag==TRUE]
    
    # get total counts unmethylated in the tumour
    UM_t <- UM_t1 + UM_t2
    #UM_t[flag==TRUE] <- round2(UM_t1[flag==TRUE]/2, digits=0) + UM_t3[flag==TRUE]
    
    # get average methyaltion rate in the tumour cells
    m_t <- (M_t)/(M_t + UM_t)
    cov_t <- M_t + UM_t
    
    # get cov_n
    cov_n <- df_mt$cov_b-cov_t
        
    # get counts methylated in the normal
    Sys.sleep(1)
    M_n <- unlist(mcmapply( function (cov, prop) rbinom(1,cov,prop) , cov=cov_n, prop=prop_n, mc.cores = n_cores)) 
    UM_n <- cov_n - M_n
    m_n <- M_n/(M_n+UM_n)
        
    # get counts methylated in the bulk
    M_b <- M_t + M_n
    UM_b <- UM_t + UM_n
    m_b <- M_b/(M_b+UM_b)
     
    # build HDI results df
    n <- nrow(df_mt) 
    HDI <- data.frame(low_mt=numeric(n), high_mt=numeric(n), 
                      low_mb=numeric(n), high_mb=numeric(n),
                      low_mn=numeric(n), high_mn=numeric(n))
    # get mt HDI
    flag_n <- !is.na(m_n) 
    flag_b <- !is.na(m_b) 
    tmp <- as.data.frame(t(mcmapply(function(M_b,UM_b,M_n,UM_n,p,CN,flag_n,flag_b) {
          x <- c(lower = as.numeric(NA), upper = as.numeric(NA))
          if(flag_n==F & flag_b == T){
            x <- HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M_b+1,shape2=UM_b+1)
          }
          if(flag_n==T){
            x <- HDIofMCMC(M_b,UM_b,M_n,UM_n,p,CN,0.99)
          }
          return(x)
    },M_b=M_b,UM_b=UM_b,M_n=M_n,UM_n=UM_n,p=p,CN=CN,flag_n=flag_n,flag_b=flag_b,mc.cores=n_cores)))
    HDI[,1:2] <- tmp
    rm(tmp,  flag_b, flag_n)
 
    # get mb HDI
    HDI[,3:4]<- as.data.frame(t(mcmapply(function(M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M+1,shape2=UM+1),mc.cores=n_cores,M=M_b,UM=UM_b)))
    
    # get mn HDI
    HDI[,5:6]<- as.data.frame(t(mcmapply(function(M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M+1,shape2=UM+1),mc.cores=n_cores,M=M_n,UM=UM_n)))
    
    # get overlaps
    ov_t <- t(mapply(function(w,x,y,z) overlaps(w,x,y,z), w=HDI[,1], x=HDI[,2], y=HDI[,5], z=HDI[,6]))
    ov_b <- t(mapply(function(w,x,y,z) overlaps(w,x,y,z), w=HDI[,2], x=HDI[,3], y=HDI[,5], z=HDI[,6]))
        
    dt <- as.data.table(HDI)
    dt$ov_t <- as.logical(ov_t) ; dt$ov_b <- as.logical(ov_b)
    rm(ov_t, ov_b)
      
    # annotate methylation rates and coverages
    dt$m_t <- m_t ; dt$cov_t <- cov_t  
    dt$m_b <- m_b ; dt$cov_b <- cov_n+cov_t 
    dt$m_n <- m_n ; dt$cov_n <- cov_n 
    
    # store CN
    dt$CN <- CN
    dt$nA <- nA
    dt$nB <- nB
    dt[ , multi := paste(nA, nB, sep="_")]
    
    # Get ground truth tumour-normal DMPs
    dt$n_gound_state <- ifelse(selector==0, "M", "UM")
    dt$t1_gound_state <- ifelse(selector1==0, "M", "UM")
    dt$t2_gound_state <- ifelse(selector2==0, "M", "UM")
    # dt$t3_gound_state <- as.character(NA)
    # dt$t3_gound_state[flag==TRUE] <- ifelse(selector3==0, "M", "UM")
    
    # Get ground truth tumour-normal DMPs
    dt$DMP_true_status <- character(length=nrow(dt))
    dt$DMP_true_status[selector==selector1 & selector==selector2] <- "no" 
    dt$DMP_true_status[selector!=selector1 | selector!=selector2] <- "yes"
    #dt$DMP_true_status[selector!=selector3 & flag == TRUE] <- "yes"
    
    # # assess percent of non-DM CpGs priors (i.e. selector=selector1)
    # # which end up with differences > 0.2 after sampling
    # sum(abs(dt$m_t-dt$m_n)<=0.2 & dt$DMP_true_status=="no", na.rm=TRUE) / sum(dt$DMP_true_status=="no")
    # 
    # # get the number of CpG with M ground state CpGs sampled as M
    # (sum(selector==0 & m_n > 0.8, na.rm=TRUE)) / sum(selector==0)
    # 
    # # get the number of CpG with UM ground state CpGs sampled as M
    # (sum(selector==1 & m_n < 0.2, na.rm=TRUE)) / sum(selector==1)
 
    # remove NAs 
    dt <- dt[!is.na(dt$m_t)&!is.na(dt$m_n)&!is.na(dt$m_b),]
    dt <- dt[dt$DMP_true_status%in%c("no", "yes"),]
    
    # # prepare object to store DMP calls
    # n<-nrow(dt)
    # result <- cbind(numeric(n), numeric(n))
    # 
    # # define variables
    # alpha_n=round2((dt$m_n*dt$cov_n), digits=0)+0.5
    # beta_n=round2(dt$cov_n-(dt$m_n*dt$cov_n), digits=0)+0.5
    # alpha_t=round2((dt$m_t*dt$cov_t), digits=0)+0.5
    # beta_t=round2(dt$cov_t-(dt$m_t*dt$cov_t), digits=0)+0.5
    # alpha_b=round2((dt$m_b*dt$cov_b), digits=0)+0.5
    # beta_b=round2(dt$cov_b-(dt$m_b*dt$cov_b), digits=0)+0.5
    # 
    # # get bulk tumour-normal DMP proabilities
    # result[,1] <- mcmapply(function(alpha_n,beta_n,alpha_b,beta_b) {
    #   prob_hyper<-NULL
    #   prob_hyper<-h(alpha_n = alpha_n, beta_n=beta_n, alpha_b=alpha_b, beta_b=beta_b)
    #   if(is.null(prob_hyper)){prob_hyper<-NA}
    #   return(prob_hyper)
    # },
    # alpha_n=alpha_n,beta_n=beta_n,alpha_b=alpha_b,beta_b=beta_b,mc.cores=n_cores)
    # 
    # # get tumour-normal DMP calls
    # dt$DMP_bn <- (result[,1]>0.99|result[,1]<0.01)&abs(dt$m_b-dt$m_n)>0.2
    # dt$DMP_tn <- (result[,1]>0.99|result[,1]<0.01)&abs(dt$m_t-dt$m_n)>0.2
    #
    # # Get TNs, TPs, FPs and FNs
    # dt$DMP_pure_status <- ifelse(dt$DMP_true_status=="yes"&dt$DMP_tn==F, "FN", 
    #                          ifelse(dt$DMP_true_status=="yes"&dt$DMP_tn==T, "TP",                                                                                         
    #                          ifelse(dt$DMP_true_status=="no"&dt$DMP_tn==F, "TN", 
    #                          ifelse(dt$DMP_true_status=="no"&dt$DMP_tn==T, "FP", NA))))
    # dt$DMP_bulk_status <- ifelse(dt$DMP_true_status=="yes"&dt$DMP_bn==F, "FN", 
    #                          ifelse(dt$DMP_true_status=="yes"&dt$DMP_bn==T, "TP",  
    #                          ifelse(dt$DMP_true_status=="no"&dt$DMP_bn==F, "TN", 
    #                          ifelse(dt$DMP_true_status=="no"&dt$DMP_bn==T, "FP", NA))))    
    
    # store output
    dt$patient <- patient_id
    dt$sample <- sample_id
    # results_list[[counter]] <- dt
      
    rm(n, HDI, CN, nA, nB)
    rm(list=ls(pattern = "(alpha)|(beta)|(sel)|(M_)|(m_)|(cov)|(prop_)"))
    
    # get multiplicitity of copy number segments outside X
    dt_top_multi <- data.table(sort(dt[, table(multi)], decreasing = T))
    colnames(dt_top_multi) <- c("multip", "counts")
    
    # keep top CNs
    dt_top_multi <- dt_top_multi[counts>1E5, ]
    x <- dt_top_multi[, as.numeric(unlist(str_split(dt_top_multi$multip, pattern = "_")))]
    for(i in 1:(length(x)%/%2)) {
      if(x[i*2] + x[(i*2)-1] > 2 | x[(i*2)-1] > 2){
        dt_top_multi[i, multip := NA] 
      }
    }
    dt_top_multi <- dt_top_multi[!is.na(multip), ]
    
    #  add results columns
    dt_top_multi[, c("b0", "expected_b0", "observed_b0", "b1", "expected_b1", "observed_b1", 
                     "t0", "expected_t0", "observed_t0", "t1", "expected_t1", "observed_t1") := 0] 
    dt_top_multi <- data.frame(dt_top_multi, stringsAsFactors = FALSE)
    top_multi <- sort(dt_top_multi$multip)
    cat("Top CpG copy number states :", sep="\n")
    cat(top_multi, sep="\n")
    
    if(length(top_multi)>0){
    # count iteration
    counter = counter + 1
    
    for(multips in top_multi){
        # get CN
        nA <- as.numeric(substr(multips, 1, 1))
        nB <- as.numeric(substr(multips, 3, 3))
        CN <- nA + nB
          
        # Get normal copy number at heterozygous SNPs
        CN_norm <- 2 
          
        # extract relevant simulated CpGs
        dt_tmp <- dt[as.character(multi)==as.character(multips),]
          
        # extract relevant observed CpGs data
        df_mt_sub <- df_mt[df_mt$multi == paste(nA, nB, sep="+"),]
        df_mt_sub0 <- df_mt_sub[df_mt_sub$m_n_high < 0.2,]
        df_mt_sub1 <- df_mt_sub[df_mt_sub$m_n_low > 0.8,]
        
        # get ratio of non-DM to potentially DM
        n_0 <- nrow(df_mt_sub0)
        n_1 <- nrow(df_mt_sub1)
   
        # get ratio of non-DM to potentially DM
        ratio_0 <- sum(df_mt_sub0$m_b > 0.25, na.rm=TRUE) / n_0
        ratio_1 <- sum(df_mt_sub1$m_b < 0.75, na.rm=TRUE) / n_1
  
        # get expected peaks in m_b
        # for hypermethylated CpGs (mn=0)
        expected_b0 =  (OC*p*CN+UC*CN_norm*(1-p))/(p*CN+CN_norm*(1-p))
          
        # for hypomethylated CpGs (mn=1)
        expected_b1 = (UC*p*CN+OC*CN_norm*(1-p))/(p*CN+CN_norm*(1-p))
  
        # get expected peaks in m_t
        # for hypermethylated CpGs (mn=0)
        expected_t0 = (expected_b0*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*UC))/(p*CN)
   
        # for hypomethylated CpGs (mn=1)
        expected_t1 = (expected_b1*(p*CN+CN_norm*(1-p))-(CN_norm*(1-p)*OC))/(p*CN)
  
        # get observed peaks in m_b
        # for bi-allelic hypermethylated CpGs (mn=0)
        x = df_mt_sub0[df_mt_sub0$m_t > 0.7  & df_mt_sub0$m_b != 1, c("m_b", "m_t")]
        
        # find mean from combined data
        m <- betareg(sort(x$m_b) ~ 1)
        observed_b0 <- unname(plogis(coef(m)[1])) # bound by [0,1]
  
        # find mean
        m <- glm(sort(x$m_t) ~ 1, family = gaussian(link = "identity"))
        observed_t0 <- unname((coef(m)[1])) # bound by [0,1]
        
        # for bi-allelic hypomethylated CpGs (mn=1)
        x = df_mt_sub1[df_mt_sub1$m_t < 0.3 & df_mt_sub1$m_b != 0, c("m_b", "m_t")]
         
        # find mean
        m <- betareg(sort(x$m_b) ~ 1)
        observed_b1 <- unname(plogis(coef(m)[1])) # bound by [0,1]
        
        # find mean
        m <- glm(sort(x$m_t) ~ 1, family = gaussian(link = "identity"))
        observed_t1 <- unname((coef(m)[1])) # bound by [0,1]
        
        # get bulk dist at ASCN = multi and for UM normal 
        n <- n_0
        n1 <- round2(n_0*ratio_0, digits = 0)
        n2 <- n1
        # n1 <- round2(n_0*ratio_0*(nA/CN), digits = 0)
        # if(as.numeric(nA)==2 & as.numeric(nB)==0){
        #   n2 <- round2(n_0*ratio_0*(1/CN), digits = 0)
        # } else {
        #   n2 <- round2(n_0*ratio_0*(nB/CN), digits = 0)
        # }
        db1 <- data.table(
            x = c(sample(dt_tmp$m_b[dt_tmp$n_gound_state=="UM" & 
                                    dt_tmp$t1_gound_state=="UM" &
                                    dt_tmp$t2_gound_state=="UM" & 
                                    dt_tmp$m_b!=1&dt_tmp$m_b!=0], n, replace=TRUE),
                  sample(dt_tmp$m_b[dt_tmp$n_gound_state=="UM" & 
                                    dt_tmp$t1_gound_state=="M" &
                                    dt_tmp$t2_gound_state=="M" &
                                    dt_tmp$m_b!=1&dt_tmp$m_b!=0], n1, replace=TRUE),
                  sample(dt_tmp$m_b[dt_tmp$n_gound_state=="UM" &
                                    dt_tmp$t1_gound_state!=dt_tmp$t2_gound_state &
                                    dt_tmp$m_b!=1&dt_tmp$m_b!=0], n2, replace=TRUE)),
          state = c(rep("ground",n), rep("dm_bi", n1), rep("dm_mono", n2))) 
          db2 <- db1
         
          ## separate histograms for both clusters
          pb1 <- ggplot(db1) + theme_minimal() +
              xlab("sim bulk tumour methylation rate")  +
            theme(legend.position = "none") +
              geom_histogram(aes(x=x, y=..count.., fill=state), 
                             binwidth = 0.02, position = "stack", alpha=0.75) + 
              geom_vline(aes(xintercept = expected_b0), size=0.5) +
              geom_vline(aes(xintercept = observed_b0), size=0.5, linetype="21") +
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80))) 
          tmp <- ggplot_build(pb1)$data[[1]]
          max_y <- 3*max(c(tmp[tmp$group==2&tmp$xmax>=0.15,"count"], 
                           tmp[tmp$group==3&tmp$xmax>=0.15,"count"]))
          suppressWarnings(pb1 <- pb1 + coord_cartesian(ylim=c(0, max_y)))
          
          # plot histogram 
          pb2 <- ggplot(db2) + theme_minimal()  +
            theme(legend.position = "none") + 
            xlab("sim bulk tumour methylation rate") +
            geom_histogram(aes(x=x, y=..count.. , fill=state), 
                           binwidth = 0.02, position = "fill", alpha=0.75) + 
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80)))
          
          # find mean from combined data
          m <- betareg(x ~ 1, data = db2[(db2$state=="ground"&db2$x>0.2)|db2$state!="ground",])
    
          # extract mu
          b0 <- unname(plogis(coef(m)[1])) # bound by [0,1]
     
          # get b1
          n <- n_1
          n1 <- round2(n_0*ratio_0, digits = 0)
          n2 <- n1
          # n1 <- round2(n_1*ratio_1*(nA/CN), digits = 0)
          # if(as.numeric(nA)==2 & as.numeric(nB)==0){
          #   n2 <- round2(n_1*ratio_1*(1/CN), digits = 0)
          # } else {
          #   n2 <- round2(n_1*ratio_1*(nB/CN), digits = 0)
          # }
  
          db3 <- data.table(
            x = c(sample(dt_tmp$m_b[dt_tmp$n_gound_state=="M" & 
                                      dt_tmp$t1_gound_state=="M" &
                                      dt_tmp$t2_gound_state=="M" & 
                                      dt_tmp$m_b!=1&dt_tmp$m_b!=0], n, replace=TRUE),
                  sample(dt_tmp$m_b[dt_tmp$n_gound_state=="M" & 
                                      dt_tmp$t1_gound_state=="UM" &
                                      dt_tmp$t2_gound_state=="UM" &
                                      dt_tmp$m_b!=1&dt_tmp$m_b!=0], n1, replace=TRUE),
                  sample(dt_tmp$m_b[dt_tmp$n_gound_state=="M" &
                                      dt_tmp$t1_gound_state!=dt_tmp$t2_gound_state &
                                      dt_tmp$m_b!=1&dt_tmp$m_b!=0], n2, replace=TRUE)),
            state = c(rep("ground",n), rep("dm_bi", n1), rep("dm_mono", n2))) 
          db4 <- db3
          
          ## separate histograms for both clusters
          pb3 <- ggplot(db3) + theme_minimal() +
            xlab("sim bulk tumour methylation rate")  +
            theme(legend.position = "none") +
            geom_histogram(aes(x=x, y=..count.., fill=state), 
                           binwidth = 0.02, position = "stack", alpha=0.75) + 
            geom_vline(aes(xintercept = expected_b1), size=0.5) +
            geom_vline(aes(xintercept = observed_b1), size=0.5, linetype="21") +
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80))) 
          tmp <- ggplot_build(pb3)$data[[1]]
          max_y <- 3*max(c(tmp[tmp$group==2&tmp$xmax<0.85,"count"], 
                           tmp[tmp$group==3&tmp$xmax<0.85,"count"]))
          suppressWarnings(pb3 <- pb3 + coord_cartesian(ylim=c(0, max_y)))
          
          # plot histogram 
          pb4 <- ggplot(db4) + theme_minimal() + 
            xlab("sim bulk tumour methylation rate") +
            geom_histogram(aes(x=x, y=..count.. , fill=state), 
                           binwidth = 0.02, position = "fill", alpha=0.75) + 
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80)))
          
          # Generate bimodal data
          m <- betareg(x ~ 1, data = db4[(db4$state=="ground"&db4$x<0.85)|db4$state!="ground",])
          
          # extract mu
          b1 <- unname(plogis(coef(m)[1])) # bound by [0,1]
          
          # make plot grid 
          gglist1 <- list(pb1, pb2)
          gglist2 <- list(pb3, pb4)
          my_layout <- cbind(1,1,1,2,2,2,3,3,3,4,4,4,4)
          gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
        
          # save plots
          f_nm <- paste0(path_patient, "Methylation/", sample_id ,
                         "/re_segmented/simulated_", nA, "_", nB,".pdf")
          ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height = 2.5, units = "in")
          
          ## color scheme observed data
          colss<-ifelse(multips%in%c("2_0","4_0"),"dodgerblue",
                 ifelse(multips%in%c("2_1","4_2"),"orange3",
                 ifelse(multips%in%c("1_1","2_2"),"firebrick", 
                 ifelse(multips=="1_0", "grey65", 
                 ifelse(multips=="3_0", "khaki", "mediumpurple3")))))
          
          ## plot observed
          pb5 <- ggplot(df_mt_sub0) + theme_minimal() +
            xlab("obs bulk tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=m_b, y=..count..,  fill=multips),  
                           binwidth = 0.02, position = "stack", alpha=0.75) + 
            geom_vline(aes(xintercept = expected_b0), size=0.5) +
            geom_vline(aes(xintercept = observed_b0), size=0.5, linetype="21") +
            scale_fill_manual(name="ASCN", values=colss)
          tmp <- ggplot_build(pb5)$data[[1]]
          max_y <- 3*max(tmp[tmp$x>=0.15,"count"])
          suppressWarnings(pb5 <- pb5 + coord_cartesian(ylim=c(0, max_y), xlim=c(0, 1)))
          suppressWarnings(pb1 <- pb1 + theme(legend.position = "bottom") +
                       coord_cartesian(ylim=c(0, max_y), xlim=c(0, 1)))
          
          ## plot observed
          pb6 <- ggplot(df_mt_sub1) + theme_minimal() +
            xlab("obs bulk tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=m_b, y=..count..,  fill=multips),  
                           binwidth = 0.02, position = "stack", alpha=0.75) + 
            geom_vline(aes(xintercept = expected_b1), size=0.5) +
            geom_vline(aes(xintercept = observed_b1), size=0.5, linetype="21") +
            scale_fill_manual(name="ASCN", values=colss)
          tmp <- ggplot_build(pb6)$data[[1]]
          max_y <- 3*max(tmp[tmp$x<=0.85,"count"])
          pb6 <- pb6 + coord_cartesian(ylim=c(0, max_y), xlim=c(0, 1))
          
          suppressWarnings(pb3 <- pb3 + theme(legend.position = "bottom")  +
                                        coord_cartesian(ylim=c(0, max_y), xlim=c(0, 1))
          )
          gglist1 <- list(pb1, pb5)
          gglist2 <- list(pb3, pb6)
          my_layout <- cbind(1,2,3,4)
          gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
          
          f_nm <- paste0(path_patient, "Methylation/", sample_id ,
                         "/re_segmented/simulated_vs_observed_", nA, "_", nB,".pdf")
          ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height = 2.5, units = "in")
          
          # resample normal from decon
          UM <- df_mt$m_n[df_mt$m_n_high < 0.2]
          M <- df_mt$m_n[df_mt$m_n_low > 0.8]
          dt_tmp$m_n_new <- ifelse(dt_tmp$n_gound_state=="UM", sample(UM, 1), sample(M, 1))
          dt_tmp$m_t_new <- (dt_tmp$m_b * (p * CN  + (1-p) * 2) - 
                             dt_tmp$m_n_new  * (1-p) * 2) / (p * CN)
    
          # get t0
          n <- n_0
          n1 <- round2(n_0*ratio_0, digits = 0)
          n2 <- n1
          dt1 <- data.table(
              x = c(sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="UM" & 
                                      dt_tmp$t1_gound_state=="UM" &
                                      dt_tmp$t2_gound_state=="UM"], n, replace=TRUE),
                    sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="UM" & 
                                      dt_tmp$t1_gound_state=="M" &
                                      dt_tmp$t2_gound_state=="M"], n1, replace=TRUE),
                    sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="UM" &
                                      dt_tmp$t1_gound_state!=
                                      dt_tmp$t2_gound_state], n2, replace=TRUE)),
            state = c(rep("ground",n), rep("dm_bi", n1), rep("dm_mono", n2))) 
     
          #Gaussian Mixed Model
          m <- glm(x ~ 1 , data = dt1[dt1$x>0.75,], family = gaussian(link = "identity"))
          t0 <- unname(coef(m)[1])
    
          ## plot observed
          pt1 <- ggplot(dt1) + theme_minimal() +
            xlab("sim pure tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=x, y=..count.., fill=state), 
                           binwidth = 0.02, position = "stack", alpha=0.75) +
            geom_vline(aes(xintercept = expected_t0), size=0.5) +
            geom_vline(aes(xintercept = observed_t0), size=0.5, linetype="21") +
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80))) 
    
          tmp <- ggplot_build(pt1)$data[[1]]
          max_y <- 3*max(tmp[tmp$x>=0.5,"count"])
          suppressWarnings(
          pt1 <- pt1 + coord_cartesian(ylim=c(0, max_y), xlim=c(-0.15, 1.5))
          )
          
          ## plot observed
          pt2 <- ggplot(df_mt_sub0) + theme_minimal() +
            xlab("sim pure tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=m_t, y=..count..,  fill=multips),
                           binwidth = 0.02, position = "stack", alpha=0.75) +
            geom_vline(aes(xintercept = expected_t0), size=0.5) +
            geom_vline(aes(xintercept = observed_t0), size=0.5, linetype="21") +
            scale_fill_manual(name="ASCN", values=colss)
    
          tmp <- ggplot_build(pt2)$data[[1]]
          max_y <- 3*max(tmp[tmp$x>=0.5,"count"])
          
          suppressWarnings(
          pt2 <- pt2 + coord_cartesian(ylim=c(0, max_y), xlim=c(-0.15, 1.5))
          )
          
          # get t1
          n <- n_1
          n1 <- round2(n_1*ratio_1, digits = 0)
          n2 <- n1
    
          dt2 <- data.table(
              x = c(sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="M" & 
                                        dt_tmp$t1_gound_state=="M" &
                                        dt_tmp$t2_gound_state=="M"], n, replace=TRUE),
                    sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="M" & 
                                        dt_tmp$t1_gound_state=="UM" &
                                        dt_tmp$t2_gound_state=="UM"], n1, replace=TRUE),
                    sample(dt_tmp$m_t_new[dt_tmp$n_gound_state=="M" &
                                        dt_tmp$t1_gound_state!=
                                        dt_tmp$t2_gound_state], n2, replace=TRUE)),
              state = c(rep("ground",n), rep("dm_bi", n1), rep("dm_mono", n2))) 
      
          #Gaussian Mixed Model
          m <- glm(x ~ 1 , data = dt2[dt2$x<0.25,], family = gaussian(link = "identity"))
          t1 <- unname(coef(m)[1])
    
          ## plot observed
          pt3 <- ggplot(dt2) + theme_minimal() +
            xlab("sim pure tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=x, y=..count.., fill=state), 
                           binwidth = 0.02, position = "stack", alpha=0.75) +
            geom_vline(aes(xintercept = expected_t1), size=0.5) +
            geom_vline(aes(xintercept = observed_t1), size=0.5, linetype="21") +
            scale_fill_manual(name="", values=c("ground"=hcl(0, 120, 80),
                                                "dm_bi"=hcl(240, 120, 80),
                                                "dm_mono"=hcl(85, 120, 80))) 
    
          tmp <- ggplot_build(pt3)$data[[1]]
          max_y <- 3*max(tmp[tmp$x<0.5,"count"])
          pt3 <- pt3 + coord_cartesian(ylim=c(0, max_y), xlim=c(-0.5, 1.15))
          
          ## plot observed
          pt4 <- ggplot(df_mt_sub1) + theme_minimal() +
            xlab("sim pure tumour methylation rate")  +
            theme(legend.position = "bottom") +
            geom_histogram(aes(x=m_t, y=..count..,  fill=multips),
                           binwidth = 0.02, position = "stack", alpha=0.75) +
            geom_vline(aes(xintercept = expected_t1), size=0.5) +
            geom_vline(aes(xintercept = observed_t1), size=0.5, linetype="21") +
            scale_fill_manual(name="ASCN", values=colss)
    
          tmp <- ggplot_build(pt4)$data[[1]]
          max_y <- 3*max(tmp[tmp$x<0.5,"count"])
          pt4 <- pt4 + coord_cartesian(ylim=c(0, max_y), xlim=c(-0.5, 1.15))
     
          gglist1 <- list(pt1, pt2)
          gglist2 <- list(pt3, pt4)
          my_layout <- cbind(1,2,3,4)
          gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
          
          f_nm <- paste0(path_patient, "Methylation/", sample_id ,
                         "/re_segmented/t_simulated_vs_observed_", nA, "_", nB,".pdf")
          ggsave(file = f_nm, plot=gd, device = "pdf", width=12, height = 2.5, units = "in")
          
          # store result
          cols <- c("b0", "expected_b0", "observed_b0", "b1", "expected_b1", "observed_b1", 
                    "t0", "expected_t0", "observed_t0", "t1", "expected_t1", "observed_t1")
          y <- c(b0, expected_b0, observed_b0, b1, expected_b1, observed_b1,
                 t0, expected_t0, observed_t0, t1, expected_t1, observed_t1)  
          dt_top_multi[dt_top_multi$multip==multips, cols] <- y 
      }
      dt_top_multi$patient <- patient_id
      dt_top_multi$sample <- sample_id
      dt_top_multi_list[[counter]] <- dt_top_multi
    }
    # clean up
    rm(list=ls(pattern = "(p[bt][1-4])|(d[bt][1-4])|([btn][012])"))
    } 
}

save(dt_top_multi_list , file = "~/rdeocd/txMeth/CAMDAC_bulk_and_pure_DM_peaks_observed_vs_expected_vs_simulated.RData")

dt_top_multi <- rbindlist(dt_top_multi_list)

tmp <- dt_top_multi
p1 <- ggplot(tmp, aes(x=observed_t1, y=expected_t1)) + 
      geom_point(shape=16, size=1, alpha=0.5, col="firebrick")+
      geom_point(aes(x=observed_t0, y=expected_t0), shape=16, size=1, alpha=0.5, col="dodgerblue")+
      xlab("Observed") + ylab("Expected") + theme_minimal() 

p2 <- ggplot(tmp, aes(x=observed_t1-expected_t1, y=..density..)) + 
      geom_histogram(alpha=0.25, fill="firebrick", binwidth=0.01) +
      geom_histogram(aes(x=observed_t0-expected_t0), alpha=0.25, fill="dodgerblue", binwidth=0.01) +
      xlab("Obs - Exp") + ylab("Density") + theme_minimal() + 
      coord_cartesian(xlim = c(-0.5, 0.5)) +
      geom_text(data=tmp[1,], aes(x = 0.15, y = 5), size=3, label=paste0("\ncor = ", 
                round2(cor(c(tmp[, observed_t1],tmp[, observed_t0]), 
                       c(tmp[, expected_t1],tmp[, expected_t0])), 
                       digits = 2)), hjust=0) 

p3 <- ggplot(tmp, aes(x=observed_t1, y=t1)) + 
  geom_point(shape=16, size=1, alpha=0.5, col="firebrick")+
  geom_point(aes(x=observed_t0, y=t0), shape=16, size=1, alpha=0.5, col="dodgerblue")+
  xlab("Observed") + ylab("Simulation") + theme_minimal() 

p4 <- ggplot(tmp, aes(x=observed_t1-t1, y=..density..)) + 
  geom_histogram(alpha=0.25, fill="firebrick", binwidth=0.01) +
  geom_histogram(aes(x=observed_t0-t0), alpha=0.25, fill="dodgerblue", binwidth=0.01) +
  xlab("Obs - Simulation") + ylab("Density") + theme_minimal() + 
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  geom_text(data=tmp[1,], aes(x = 0.15, y = 5), size=3, 
            label=paste0("\ncor = ", round2(cor(c(tmp[, observed_t1],tmp[, observed_t0]), 
                                                c(tmp[, t1],tmp[, t0])), 
                                            digits = 2)), 
            hjust=0) 

p5 <- ggplot(tmp, aes(x=t1, y=expected_t1)) + 
  geom_point(shape=16, size=1, alpha=0.5, col="firebrick")+
  geom_point(aes(x=t0, y=expected_t0), shape=16, size=1, alpha=0.5, col="dodgerblue")+
  xlab("Simulation") + ylab("Expected") + theme_minimal() 

p6 <- ggplot(tmp, aes(x=t1-expected_t1, y=..density..)) + 
  geom_histogram(alpha=0.25, fill="firebrick", binwidth=0.01) +
  geom_histogram(aes(x=t0-expected_t0), alpha=0.25, fill="dodgerblue", binwidth=0.01) +
  xlab("Simulation - Exp") + ylab("Density") + theme_minimal() + 
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  geom_text(data=tmp[1,], aes(x = 0.15, y = 5), size=3, 
            label=paste0("\ncor = ", round2(cor(c(tmp[, t1],tmp[, t0]), 
                                                c(tmp[, expected_t1],tmp[, expected_t0])), 
                                            digits = 2)), 
            hjust=0) 

# make grid 
gglist1 <- list(p1, p2, p3) ; gglist2 <- list(p4, p5, p6) 
my_layout <- cbind(1, 2, 3, 4, 5, 6)
gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)

# save grid
f_nm <- "~/rdeocd/txMeth_analysis/CAMDAC_pure_DM_peaks_observed_vs_expected_vs_simulated.pdf"
ggsave(file = f_nm, plot=gd, device = "pdf", width=14, height=2.5, units = "in")
rm(gglist1, gglist2, my_layout, gd, f_nm, p1, p2, p3, p4, p5, p6)

# load purity info
df_purity <- data.table(read_tsv("~/rdeocd/txMeth/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"))
df_purity <- df_purity[!grepl("LN",Sample),]

# set patient ids
tmp$LTX_ID <- NULL
tmp$LTX_ID <- paste(tmp$patient,tmp$sample,sep = "_")

# assign rho to dt
tmp$rho <- df_purity$rho_reseg[match(tmp$LTX_ID, df_purity$LTX_ID)]

# order by purity 
tmp2 <- tmp[, .(mean_rho = mean(rho, na.rm=TRUE)), by = .(LTX_ID)]
tmp2[, order := order(mean_rho)]


# set sample ids as ordered factor
tmp$LTX_ID <- factor(tmp$LTX_ID, levels = tmp2$LTX_ID[tmp2$order], 
                     ordered = T)  

# plot rho & rho_reseg per patient
ggplot(tmp, aes(x=LTX_ID, y=(observed_t1-expected_t1)))+ theme_classic() +
  geom_violin(trim = FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', position = "dodge", dotsize = .85, binwidth = .005, alpha=.5) + 
  # scale_colour_manual(name="Method", values = c("CAMDAC\nmulti"="darkorange","CAMDAC\nsingle"="firebrick")) +
  # scale_fill_manual(name="Method", values = c("CAMDAC\nmulti"="darkorange","CAMDAC\nsingle"="firebrick")) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) + 
  scale_y_continuous(name="Obs - Exp",limits = c(-0.3,0.3), breaks = seq(-0.3,0.3, 0.1), 
                     labels = c("-0.3", "-0.2", "-0.1", "0.0", "0.1", "0.2", "0.3")) +
  ggsave(filename = "~/rdeocd/txMeth_analysis/CAMDAC_pure_DM_peaks_obs_vs_exp_per_sample.pdf", 
         device = "pdf", units = "in", width = 12, height = 3)
