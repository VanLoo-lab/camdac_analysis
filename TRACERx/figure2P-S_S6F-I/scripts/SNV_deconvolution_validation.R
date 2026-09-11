#' @description CAMDAC deconvolution validation with SNVs in regions 
#' of tumour loss of heterozygosity (LOH). 
#' 
#' We also validate the adjacent normal by comparing to the normal 
#' infiltrates methylation taking the wild type reads in regions
#' of LOH. This could be extended to SNPs in regions of LOH as all
#' reads from the lost allele would also exclusively come from normal
#' cells.
#' 
#' 
#' @author Lili Larose Cadieux

# Load library
x <- c("data.table", "stringr", "ggplot2", "gridExtra", "scales",
       "ggpubr", "matrixStats", "parallel")

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

extend <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

lm_eqn <- function(df, x, y) {
  m <- lm(y ~ 0 + x, df);
  eq <- substitute(italic(y) == a %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
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
intervalWidth <- function(lowTailPr, ICDFname, credMass, ... ) {
  ICDFname(credMass + lowTailPr, ... ) - ICDFname(lowTailPr, ... )
}

HDIofICDF <- function(ICDFname, credMass = 0.99 , tol = 1e-4, ... ) {

  incredMass <- 1.0 - credMass

  optInfo <-
    stats::optimize(
      f = intervalWidth,
      interval = c(0,incredMass),
      ICDFname = ICDFname,
      credMass = credMass,
      tol = tol,
    ... )

  HDIlowTailPr <- optInfo$minimum
  vec <-
    setNames(
      object = ICDFname(c(HDIlowTailPr, credMass+HDIlowTailPr), ... ),
      nm = c("low", "high")
    )

  return(vec)
}

################################################################################
# Set paths and load inputs
################################################################################

# set parent directory
tx_dir <- "~/rdeocd/txMeth"

# set the working directory
work_dir <- "~/data/CAMDAC_paper/figure2P-S_S6F-I"
setwd(work_dir)

# set the output directory
out_dir <- file.path(work_dir, "outputs", Sys.Date())
if (!dir.exists(dirname(out_dir))) dir.create(dirname(out_dir))
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# Load TRACERx and CRC data
################################################################################

# load the TRACERx cohort phased data
f_nm <- "inputs/TRACERx_SNV_DMPs.RData"
load(f_nm)
dt_tx <- as.data.table(all_df_phased_methylation)
dt_tx$cohort <- "tracerx"
rm(all_df_phased_methylation)

# load the CRC cohort RRBS phased data
f_nm <- "inputs/CRC_2PQ_RRBS_new.rds"
tmp <- readRDS(f_nm)
dt_crc <- as.data.table(tmp)
dt_crc$cohort <- "crc_rrbs"
rm(tmp)

# load the CRC cohort WGBS phased data
#f_nm <- "inputs/CRC_2N_WGBS.RDS"
f_nm <- "inputs/CRC_2N_WGBS_VAF0.1_ALT3.RDS"
tmp <- readRDS(f_nm)
dt_crc_wgbs <- as.data.table(tmp)
dt_crc_wgbs$cohort <- "crc_wgbs"
rm(tmp)

# fill in missing info
# For clonal SNVs on all tumour copied in regions on LOH, m_2 == m_t_2
dt_crc_wgbs[, "m_t_2" := m_2]
dt_crc_wgbs[, "m_diff_tn2" := m_n - m_t_2]

################################################################################
# Annotated tumour purity
################################################################################

# load tracerx patient information table with annotated tumour purity
patient_f_nm_tx <- paste0("tx100clinicalinfo_RRBS_published_status_",
    "immune_clusters_purity_ploidy_CpG_cov.txt")
patient_info_tx <- fread(file.path(tx_dir, patient_f_nm_tx))
patient_info_tx <- patient_info_tx[!grepl("LN", Sample)]

# add purity values
dt_tx$rho <- patient_info_tx[match(dt_tx$LTX_ID, LTX_ID), "rho_reseg"]

# load CRC RRBS tumour purity
purity_dir <- "../figure1H-K/inputs/"
patient_f_nm_crc <- "CRC_RRBS_WGS_overlap.txt"
patient_info_crc <-
  fread(file.path(purity_dir, patient_f_nm_crc), sep = ",")

# add purity values
dt_crc$rho <- patient_info_crc[match(dt_crc$PT_ID, sample_id), "purity_b"]
dt_crc_wgbs$rho <-
  patient_info_crc[match(dt_crc_wgbs$PT_ID, sample_id), "purity_a"]

################################################################################
# Annotated tumour ploidy
################################################################################

# add purity values
dt_tx$ploidy <- patient_info_tx[match(dt_tx$LTX_ID, LTX_ID), "psi_reseg"]

# add purity values
dt_crc$ploidy <- patient_info_crc[match(dt_crc$PT_ID, sample_id), "ploidy_b"]
dt_crc_wgbs$ploidy <-
  patient_info_crc[match(dt_crc_wgbs$PT_ID, sample_id), "ploidy_a"]

################################################################################
# Subset columns and compute necessary variables
################################################################################

# columns of interest
cols <- c("chrom", "start", "POS_genotype", "width",
          "m_n", "M_n", "UM_n", "m_b", "m_t",
          "m_1", "m_2", "M_2", "UM_2",
          "m_t_2", "m_diff_tn2",
          "nA", "nB",
          "SNV_mcn_methylation", "SNV_mcn",
          "SNV_VAF_WGS", "SNV_VAF_methylation",
          "rho", "ploidy", "cohort", "sequencing_platform")

# combine the data
dt_tx[, "sequencing_platform" := "RRBS"]
dt_crc[, "sequencing_platform" := "RRBS"]
dt_crc_wgbs[, "sequencing_platform" := "WGBS"]
dt_all <- rbindlist(list(dt_tx[, ..cols], dt_crc[, ..cols], dt_crc_wgbs[, ..cols]))

# add m_t_2_corr column
dt_all[, "m_t_2_corr" := ifelse(m_t_2 > 1, 1, ifelse(m_t_2 < 0, 0, m_t_2))]

# Calculate HDI for the methylation rate
n <- nrow(dt_all)
M_n <- dt_all$M_n
UM_n <- dt_all$UM_n
vec <- cbind.data.frame(low = numeric(length = n), high = numeric(length = n))
vec[, 1:2] <- do.call(rbind,
  mclapply(1:n, function(i,M,UM)
          HDIofICDF(
            ICDFname = qbeta,
            credMass = 0.99,
            shape1 = M[i] + 1,
            shape2 = UM[i] + 1),
            M = M_n,
            UM = UM_n,
            mc.cores = 4))
dt_all[, c("m_n_low", "m_n_high") := vec]

# Calculate HDI for the mutant allele methylation rate
M_2 <- dt_all$M_2
UM_2 <- dt_all$UM_2
vec <- cbind.data.frame(low = numeric(length = n), high = numeric(length = n))
vec[, 1:2] <- do.call(rbind,
    mclapply(1:n, function(i,M,UM)
          HDIofICDF(
            ICDFname = qbeta,
            credMass = 0.99,
            shape1 = M[i] + 1,
            shape2 = UM[i] + 1),
            M = M_2,
            UM = UM_2,
            mc.cores = 4))
dt_all[, c("m_t_2_low", "m_t_2_high") := vec]

# note VAF difference - RRBS vs WGS
dt_all[, vaf_diff := abs(SNV_VAF_methylation - SNV_VAF_WGS)]

# remove NAs and subset to LOH regions
rows_to_keep <- !is.na(dt_all$m_t) & !is.na(dt_all$m_t_2) & dt_all$nB == 0
cat("Number of CpGs removed due to NAs or non-zero nB:", sum(!rows_to_keep), "\n")
dt_all[, table(!rows_to_keep, cohort)]
dt_trimmed <- dt_all[rows_to_keep]

# remove subclonal CpGs
#dt <- dt[SNV_mcn_methylation >= 0.5 & abs(SNV_mcn_methylation - nA) <= 0.5]
rows_to_keep <- dt_trimmed$SNV_mcn >= 0.75 & abs(dt_trimmed$SNV_mcn - dt_trimmed$nA) <= 0.25
cat("Number of CpGs removed due to subclonal SNVs:", sum(!rows_to_keep), "\n")
dt_trimmed[, table(!rows_to_keep, cohort)]
dt_trimmed <- dt_trimmed[rows_to_keep]

# remove low VAF SNVs
rows_to_keep <- dt_trimmed$SNV_VAF_WGS >= 0.1
cat("Number of CpGs removed due to low VAF:", sum(!rows_to_keep), "\n")
dt_trimmed[, table(!rows_to_keep, cohort)]
dt_trimmed <- dt_trimmed[rows_to_keep]

# remove CpGs with < 3 alt allele phased reads
rows_to_keep <- dt_trimmed$M_2 + dt_trimmed$UM_2 >= 3
cat("Number of CpGs removed due to low alt allele phased reads:", sum(!rows_to_keep), "\n")
dt_trimmed[, table(!rows_to_keep, cohort)]
dt_trimmed <- dt_trimmed[rows_to_keep]

# Check for duplicates - there should be none
if (any(duplicated(dt_trimmed))) {
  stop("Duplicated CpGs found in the data")
}

# note VAF difference - RRBS vs WGS
summary_stats <- dt_trimmed[!is.na(vaf_diff), 
        .(variable = c("min", "q25", "q50", "mean", "q75", "max"), 
        val = summary(vaf_diff)), by = .(cohort)]
summary_stats_wide <- dcast(summary_stats, cohort ~ variable, value.var = "val")
fwrite(summary_stats_wide,
  file.path(out_dir, "summary_stats_wide_VAF_diff.tsv"), sep = "\t")

# get capped methylation rates
dt_trimmed[, m_t_corr := ifelse(m_t < 0, 0, ifelse(m_t > 1, 1, m_t))]
dt_trimmed[, m_t_2_corr := ifelse(m_t_2 < 0, 0, ifelse(m_t_2 > 1, 1, m_t_2))]

# extract WT methylation rate
dt_trimmed[, m_wt := m_1]
# add immune category # LEGACY
# dt[, "immune_category" := patient_info[match(dt$LTX_ID, LTX_ID), immune_cluster]]

# loop over sequencing platforms
for (sequencing_platform in c("RRBS", "WGBS", "BS")) {

  if (sequencing_platform == "RRBS") {
    dt <- dt_trimmed[sequencing_platform == "RRBS"]
  } else if (sequencing_platform == "WGBS") {
    dt <- dt_trimmed[sequencing_platform == "WGBS"]
  } else {
    dt <- dt_trimmed
  }

  # note down the number of phased clonal SNVs in regions of LOH
  dt_cpgs_per_snvs <- dt[, .(num_CpGs = .N), by = .(chrom, POS_genotype, cohort)]
  tb1 <- dt_cpgs_per_snvs[, .(mean_num_CpGs = mean(num_CpGs),
                      total_num_CpGs = sum(num_CpGs)), by = .(cohort)]
  tb2 <- dt_cpgs_per_snvs[, .(num_snvs = .N), by = .(cohort)]

  # save the data
  f_nm <- paste0("cpgs_per_snv_", sequencing_platform, ".txt")
  fwrite(tb1, file.path(out_dir, f_nm), sep = "\t")

  # save the data
  f_nm <- paste0("snv_phased_to_cpgs_", sequencing_platform, ".txt")
  fwrite(tb2, file.path(out_dir, f_nm), sep = "\t")

  ################################################################################
  # Methylation-derived VAF validation
  ################################################################################

  # NOTE: SNV_VAF = (SNV_mcn * rho) / (nA * rho + 2 * (1 - rho))

  # compare observed vs expected epimutation copy number
  p1 <- ggplot(dt, aes(x = SNV_mcn, y = SNV_mcn_methylation)) +
    theme_classic() + xlab("Expected (WGS/WES)") +
    ylab("Observed (methylation)") +
    geom_point(size = 1, shape = 16, alpha = 0.5, col = "green3") +
    ylim(0, 6) + xlim(0, 6) +
    geom_text(data = dt[1, ],
      aes(x = 6 + 0.1, y = 6 - 0.1), size = 2,
          label = paste0("n = ", nrow(dt),
                        "\ncor = ", round2(cor(dt[, SNV_mcn],
                                                dt[, SNV_mcn_methylation]),
                        digits = 2)), hjust = 0)

  f_nm <- paste0("SNV_mcn_", sequencing_platform, "_vs_GS_points.pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p1, device = "pdf", units = "in", width = 3.75, height = 3.75)

  # compare observed vs expected epimutation copy number
  tmp <- dt[, .(SNV_mcn = SNV_mcn,
                SNV_mcn_methylation_norm = SNV_mcn_methylation / ploidy,
                ploidy = ploidy
                )]
  p1 <- ggplot(tmp, aes(x = SNV_mcn, y = SNV_mcn_methylation_norm)) +
    theme_classic() + xlab("Expected (WGS/WES)") +
    ylab("Observed (methylation)") +
    geom_point(size = 1, shape = 16, alpha = 0.5, col = "green3") +
    ylim(0, 3.5) + xlim(0, 6.5) +
    annotate("text", x = 0.5, y = 2.9, size = 2, hjust = 0,
      label = paste0("n = ", nrow(dt),
                    "\ncor = ", round2(cor(tmp[, SNV_mcn],
                                           tmp[, SNV_mcn_methylation_norm]),
                    digits = 2)))

  f_nm <- paste0("SNV_mcn_", sequencing_platform, "_vs_GS_points_normalised.pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p1, device = "pdf", units = "in", width = 3.75, height = 3.75)

  # compare observed vs expected vaf - BS vs GS
  p1 <- ggplot(dt, aes(x=SNV_VAF_WGS, y=SNV_VAF_methylation)) +
    theme_classic() + xlab("WGS/WES VAF") +
    ylab(paste0(sequencing_platform, " VAF")) +
    #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
    geom_point(size=1, shape=16, alpha=0.25, col="black", fill="black") +
    #geom_smooth(method="lm", col="grey40") +
    geom_abline(slope=1, intercept = 0, linewidth=0.5, linetype="dashed", col="grey40") +
    ylim(0,1) + xlim(0,1) +
    #scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") +
    geom_text(data=dt[1,], aes(x = 0.1, y = 0.8), size=2, 
              label=paste0("n = ", nrow(dt), "\ncor = ", round2(cor(
                    dt[, SNV_VAF_WGS], dt[, SNV_VAF_methylation]),
              digits = 2)), hjust=0)

  f_nm <- paste0("SNV_vaf_", sequencing_platform, "_vs_GS_points.pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

  # create a new data frame with only mutation copy number and VAF differences
  df <- data.frame(y = c(dt[,(SNV_mcn_methylation - SNV_mcn)]),
                  y2 = c(dt[,(SNV_VAF_methylation - SNV_VAF_WGS)]),
                  x = factor(c(rep(0, nrow(dt))), levels = c(0), ordered = TRUE), 
                  flag = c(rep("BS-GS", nrow(dt))))

  # plot mutation copy number differences as boxplot
  p1 <- ggplot(df, aes(x=x, y=y, col=flag, fill=flag)) +
      theme_classic() +
      ylab("Observed - Expected")  +
      geom_boxplot(width=.75, alpha=0.1) +
      scale_colour_manual(name="", values=c("BS-GS"="darkorange")) +
      scale_fill_manual(name="", values=c("BS-GS"="darkorange")) + 
      #stat_compare_means(comparisons = my_comparisons, size = 2) + 
      scale_x_discrete(name="Platform", breaks = 0, labels = c("BS-GS"))
  f_nm <- paste0("SNV_mcn_BS_vs_GS_box_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p1, device = "pdf", units = "in", width = 3.5, height = 2.8)

  # plot mutation VAF differences as boxplot
  p1 <- ggplot(df, aes(x=x, y=y2, col=flag, fill=flag)) +  theme_classic() +
    ylab("vaf (observed - expected)")  +
    geom_boxplot(width=.75, alpha=0.1) +
    scale_colour_manual(name="", values=c("BS-GS"="darkorange")) +
    scale_fill_manual(name="", values=c("BS-GS"="darkorange")) + 
    #stat_compare_means(comparisons = my_comparisons, size = 2) 
    scale_x_discrete(name="Platform", breaks = 0, labels = c("BS-GS"))
  f_nm <- paste0("SNV_vaf_BS_vs_GS_box_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

  # annotate SNVs at CCGG
  dt[, "dist_SNV" := abs(as.numeric(POS_genotype) - start)]
  dt[, "flag" := dist_SNV < 4 & width == 4]
  tmp <- dt[, .(CCGG = sum(flag) > 0), by = .(chrom, POS_genotype)]
  values <- tmp[match(paste0(dt$chrom, dt$POS_genotype),
                      paste0(chrom, POS_genotype)), CCGG]
  dt[, "flag2" := ifelse(is.na(values), FALSE, values)]
  df <- data.frame(y2 = c(dt[flag2 == FALSE, (SNV_VAF_methylation - SNV_VAF_WGS)],
                          dt[flag2 == TRUE, (SNV_VAF_methylation - SNV_VAF_WGS)]),
                  x = factor(c(rep(0, sum(dt$flag2 == FALSE, na.rm = TRUE)),
                                rep(1, sum(dt$flag2 == TRUE, na.rm = TRUE))),
                              levels = c(0,1), ordered = TRUE),
                  flag = c(rep("non-CCGG", sum(dt$flag2 == FALSE, na.rm = TRUE)),
                            rep("CCGG", sum(dt$flag2 == TRUE, na.rm = TRUE))))

  p1 <- ggplot(df, aes(x=x, y=y2, col=flag, fill=flag)) +  theme_classic() +
    ylab("vaf (observed - expected)")  +
    geom_boxplot(width=.75, alpha=0.1) +
    scale_colour_manual(name="", values=c("CCGG"="mediumpurple2","non-CCGG" = "grey50")) +
    scale_fill_manual(name="", values=c("CCGG"="mediumpurple2","non-CCGG" = "grey50")) + 
    #stat_compare_means(size = 2) +
    scale_x_discrete(name="Platform", breaks = c(0,1), labels = c("non-CCGG","CCGG"))
  f_nm <- paste0("SNV_vaf_BS_vs_GS_with_CCGG_context_box_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

  # compare observed vs expected vaf
  p1 <- ggplot(dt, aes(x=SNV_VAF_WGS, y=SNV_VAF_methylation)) +
    theme_classic() + xlab("WGS/WES VAF") +
    ylab(paste0(sequencing_platform, " VAF")) +
    #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
    geom_point(data=dt[flag2==FALSE], 
              aes(col=flag2, fill=flag2), size=1, shape=16, alpha=0.2) +
    geom_point(data=dt[flag2==TRUE], 
              aes(col=flag2, fill=flag2), size=1.5, shape=16, alpha=1) +
    #geom_smooth(method="lm", col="grey40") +
    scale_colour_manual(name="", values=c("FALSE"="grey50","TRUE"= "mediumpurple2"),
                        labels=c("FALSE"="non-CCGG", "TRUE"="CCGG")) +
    scale_fill_manual(name="", values=c("FALSE"="grey50","TRUE"="mediumpurple2"),
                      labels=c("FALSE"="non-CCGG", "TRUE"="CCGG")) +
    geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="grey40") +
    ylim(0,1) + xlim(0,1) +
    geom_text(data=dt[1,], aes(x = 0.1, y = 0.8), size=2, 
              label=paste0("n = ", nrow(dt), 
                          "\ncor = ", round2(cor(
                dt[, SNV_VAF_WGS], 
                dt[, SNV_VAF_methylation]), 
                digits = 2),
                "\ncor without CCGG = ", round2(cor(
                  dt[flag2==FALSE, SNV_VAF_WGS], 
                  dt[flag2==FALSE, SNV_VAF_methylation]), 
                digits = 2),
                "\ncor at CCGG = ", round2(cor(
                  dt[flag2==TRUE, SNV_VAF_WGS], 
                  dt[flag2==TRUE, SNV_VAF_methylation]), 
                digits = 2)), hjust=0)
  f_nm <-
    paste0("SNV_vaf_BS_vs_GS_points_with_CCGG_context_",
           sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.75, height = 2.8)

  ################################################################################
  # CAMDAC validation
  ################################################################################

  # set plot boundaries
  max_plus <- max(na.omit(c(dt$m_t, dt$m_t_2)))
  max_minus <- max(na.omit(c(-dt$m_t, -dt$m_t_2)))
  
  df <- data.frame(y = c(dt[,(m_t_2)], dt[,(m_t)]),
                  x = factor(c(rep(0, nrow(dt)),rep(1, nrow(dt))),
                             levels = c(0,1), ordered = TRUE), 
                  flag = c(rep("SNV", nrow(dt)),rep("CAMDAC", nrow(dt))))
  my_comparisons <- list(factor(c(0, 1), levels = c(0,1), ordered = TRUE))

  # compare the distribution of SNV- and CAMDAC-deconvolved methylation rates 
  p1 <- ggplot(df, aes(x=x, y=y, col=flag, fill=flag)) +  theme_classic() +  
    ylab("pure tumour methylation rate")  +
    geom_boxplot(width=.75, alpha=0.1) +
    scale_colour_manual(name="", values=c("SNV"="grey30", "CAMDAC"="darkorange")) +
    scale_fill_manual(name="", values=c("SNV"="grey30", "CAMDAC"="darkorange")) + 
    #stat_compare_means(comparisons = my_comparisons, size = 2) + 
    scale_x_discrete(name="Platform", breaks = seq(0, 1, 1), labels = c("SNV","CAMDAC"))
  f_nm <- paste0("SNV_decon_vs_CAMDAC_box_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p1, device = "pdf", units = "in", width = 3.5, height =2.8)

  # WGBS only ouputs capped methylation rates, this only applies to RRBS
  if (sequencing_platform == "RRBS") {
    # compare the point estimates of SNV- and CAMDAC-deconvolved methylation rates
    p1 <- ggplot(dt, aes(x=m_t_2, y=m_t)) + 
      geom_point(size=0.5, shape=16, alpha=0.25, col="navyblue") + theme_classic() + 
      geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
      #geom_smooth(method="lm", formula = y ~ 0 + x) +
      ylim(-max_minus,max_plus) + xlim(0,1) + 
      geom_text(data=dt[1,], aes(x = 0, y = 1.35), size=2, 
                label=paste0("n = ", nrow(dt), 
                            "\ncor = ", round2(cor(dt[, m_t], 
                                          dt[, m_t_2]), digits = 2)), hjust=0)
    f_nm <- paste0("SNV_decon_vs_CAMDAC_decon_", sequencing_platform, ".pdf")
    ggsave(filename = file.path(out_dir, f_nm),
          plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

    # evalute the distribution of differences
    p1 <- ggplot(dt, aes(x=m_t_2-m_t, y=..density..)) +
      geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() +
      geom_text(data=dt[1,], aes(x = -max_minus+0.5, y = 3), size=2,
                label=paste0(
                  "n = ", nrow(dt),
                  "\ncor = ",
                    round2(
                      cor(dt[, m_t], dt[, m_t_2]),
                      digits = 2)
                ), hjust=0)
    f_nm <- paste0("SNV_decon_vs_CAMDAC_decon_density_",
      sequencing_platform, ".pdf")
    ggsave(filename = file.path(out_dir, f_nm),
          plot=p1, device = "pdf", units = "in", width = 3, height = 2.8)
  }

  # compare the point estimates of SNV- and bulk tumour methylation rates
  # colour by |m_diff_tn2| > 0.2 (red) vs <= 0.2 (grey); annotate full (black) and DM subset (red in brackets)
  tmp <- dt[abs(m_diff_tn2) > 0.2]
  dt[, dm_tn2 := abs(m_diff_tn2) > 0.2]
  p <- ggplot(dt, aes(x = m_2, y = m_b, color = dm_tn2)) + theme_classic() +
    geom_point(size = 1, shape = 16, alpha = 0.3) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                       labels = c("FALSE" = "|m_snv-m_b| <= 0.2",
                                  "TRUE" = "|m_snv-m_b| > 0.2"),
                       name = NULL) +
    ylim(0, 1.1) + xlim(0, 1) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dashed", col = "black") +
    annotate("text", x = 0, y = 1.02, size = 2, hjust = 0, color = "black",
      label = paste0("n = ", nrow(dt), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("n = ", nrow(dt), " ")), y = 1.02, size = 2, hjust = 0, color = "red",
      label = paste0("(", nrow(tmp), ")")) +
    annotate("text", x = 0, y = 0.97, size = 2, hjust = 0, color = "black",
      label = paste0("cor = ", round2(cor(dt[, m_b], dt[, m_2]), digits = 2), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("cor = ", round2(cor(dt[, m_b], dt[, m_2]), digits = 2), " ")), y = 0.97, size = 2, hjust = 0, color = "red",
      label = paste0("(", round2(cor(tmp[, m_b], tmp[, m_2]), digits = 2), ")")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.text = element_text(color = "black", size = 6))
  f_nm <- paste0("SNV_decon_vs_bulk_tn2_diff_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p, device = "pdf", units = "in", width = 3.4 + 1, height = 2.8)

  # compare the point estimates of SNV- and CAMDAC-purified tumour methylation rates
  # colour by |m_diff_tn2| > 0.2 (red) vs <= 0.2 (grey); annotate full (black) and DM subset (red in brackets)
  p1 <- ggplot(dt, aes(x = m_2, y = m_t_corr, color = dm_tn2)) + theme_classic() +
    geom_point(size = 1, shape = 16, alpha = 0.3) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                       labels = c("FALSE" = "|m_snv-m_b| <= 0.2",
                                  "TRUE" = "|m_snv-m_b| > 0.2"),
                       name = NULL) +
    ylim(0, 1.1) + xlim(0, 1) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dashed", col = "black") +
    annotate("text", x = 0, y = 1.02, size = 2, hjust = 0, color = "black",
      label = paste0("n = ", nrow(dt), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("n = ", nrow(dt), " ")), y = 1.02, size = 2, hjust = 0, color = "red",
      label = paste0("(", nrow(tmp), ")")) +
    annotate("text", x = 0, y = 0.97, size = 2, hjust = 0, color = "black",
      label = paste0("cor = ", round2(cor(dt[, m_t_corr], dt[, m_2]), digits = 2), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("cor = ", round2(cor(dt[, m_t_corr], dt[, m_2]), digits = 2), " ")), y = 0.97, size = 2, hjust = 0, color = "red",
      label = paste0("(", round2(cor(tmp[, m_t_corr], tmp[, m_2]), digits = 2), ")")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.text = element_text(color = "black", size = 6))
  f_nm <- paste0("SNV_decon_vs_CAMDAC_decon_capped_tn2_diff_",
    sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p1, device = "pdf", units = "in", width = 3.4 + 1, height = 2.8)

  # repeat comparison of the point estimates of SNV-purified and bulk tumour
  # methylation rates where the variant and adjacent normal at significant DMPs
  tmp <- dt[abs(m_diff_tn2) > 0.2 & (m_t_2_low > m_n_high | m_t_2_high < m_n_low)]
  p1 <- ggplot(tmp, aes(x=m_2, y=m_b)) + theme_classic() +
    geom_point(size=1, shape=16, alpha=0.25, col="navyblue") +
    #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
    ylim(0,1) + xlim(0,1) +
    #scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") + 
    geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
    geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2,
              label=paste0("n = ", nrow(tmp), "\ncor = ",
                          round2(cor(tmp[, m_b], tmp[, m_2]),  digits = 2)), hjust=0)
  f_nm <- paste0("SNV_decon_vs_bulk_capped_DMP_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.4, height = 2.8)

  # repeat comparison of the point estimates of SNV- and CAMDAC-purified tumour 
  # methylation rates where the variant and adjacent normal at significant DMPs     
  p1 <- ggplot(tmp, aes(x=m_2, y=m_t_corr)) + theme_classic() + 
    geom_point(size=1, shape=16, alpha=0.25, col="navyblue") +
    #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
    ylim(0,1) + xlim(0,1) + 
    #scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") + 
    geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
    geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2,
              label=paste0("n = ", nrow(tmp), "\ncor = ",
                          round2(cor(tmp[, m_t_corr], tmp[, m_2]),  digits = 2)), hjust=0)
  f_nm <- paste0("SNV_decon_vs_CAMDAC_decon_capped_DMP_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3.4, height = 2.8)

  p1 <- ggplot(tmp, aes(x=m_t_2_corr-m_t_corr, y=..density..)) +
    geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() +
    geom_text(data=dt[1,], aes(x = -max_minus+0.5, y = 3), size=2, 
              label=paste0("n = ", nrow(tmp),  "\ncor = ", 
                          round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]), digits = 2)), hjust=0)
  f_nm <- paste0("SNV_decon_vs_CAMDAC_decon_density_capped_DMP_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

  ################################################################################
  # Adjacent normal validation
  ################################################################################

  # compare normal adjacent and WT - colour by |m_diff_tn2| > 0.2 (red) vs <= 0.2 (grey); annotate full and DM subset
  dt[, m_wt := m_1]
  dt_wt <- dt[!is.na(m_wt)]
  tmp_wt <- dt_wt[abs(m_diff_tn2) > 0.2]
  p <- ggplot(dt_wt, aes(x = m_wt, y = m_n, color = dm_tn2)) + theme_classic() +
    geom_point(size = 1, shape = 16, alpha = 0.3) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"),
                       labels = c("FALSE" = "|Δm| ≤ 0.2", "TRUE" = "|Δm| > 0.2"),
                       name = NULL) +
    ylim(0, 1.1) + xlim(0, 1) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.5, linetype = "dashed", col = "black") +
    annotate("text", x = 0, y = 1.02, size = 2, hjust = 0, color = "black",
      label = paste0("n = ", nrow(dt_wt), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("n = ", nrow(dt_wt), " ")), y = 1.02, size = 2, hjust = 0, color = "red",
      label = paste0("(", nrow(tmp_wt), ")")) +
    annotate("text", x = 0, y = 0.97, size = 2, hjust = 0, color = "black",
      label = paste0("cor = ", round2(cor(dt_wt[, m_wt], dt_wt[, m_n]), digits = 2), " ")) +
    annotate("text", x = 0.04 + 0.014 * nchar(paste0("cor = ", round2(cor(dt_wt[, m_wt], dt_wt[, m_n]), digits = 2), " ")), y = 0.97, size = 2, hjust = 0, color = "red",
      label = paste0("(", round2(cor(tmp_wt[, m_wt], tmp_wt[, m_n]), digits = 2), ")")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.text = element_text(color = "black", size = 6))
  f_nm <- paste0("WT_decon_vs_adjacent_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm),
        plot = p, device = "pdf", units = "in", width = 3.4 + 1, height = 2.8)

  # explanation for lower number of CpGs in the plot below than in the SNV_mt plot
  dt[abs(m_diff_tn2) > 0.2 & (m_t_2_low > m_n_high | m_t_2_high < m_n_low)][, 
    table(!is.na(m_wt))
    ]

  # Is the different exacerbated at DMPs
  tmp <- dt_wt[abs(m_diff_tn2) > 0.2 & (m_t_2_low > m_n_high | m_t_2_high < m_n_low)]
  p <- ggplot(tmp, aes(x=m_wt, y=m_n)) + theme_classic() + 
    geom_point(size=1, shape=16, alpha=0.25, col="black") +
    ylim(0,1) + xlim(0,1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="black") +
    geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
              label=paste0("n = ", nrow(tmp), "\ncor = ", 
                          round2(cor(tmp[, m_n], tmp[, m_wt]),  digits = 2)), hjust=0)
  f_nm <- paste0("WT_decon_vs_adjacent_DMP_", sequencing_platform, ".pdf")
  ggsave(filename = file.path(out_dir, f_nm), 
        plot = p, device = "pdf", units = "in", width = 3.4, height = 2.8)

  ## Are we less accurate in immune hot tumours?
  #p1 <- ggplot(dt_wt, aes(x = m_wt - m_n, y = ..count..)) +
  #  theme_classic() + 
  #  geom_histogram(aes(col = immune_category, fill = immune_category), alpha=0.25, bins = 50) +
  #  xlim(-1,1) 
  #ggsave(filename = file.path(out_dir, "WT_decon_vs_adjacent_by_immune_cluster.pdf"), 
  #       plot = p1, device = "pdf", units = "in", width = 4.5, height = 2.8)
  ## Legacy, no observable difference in TRACERx.
}

# END