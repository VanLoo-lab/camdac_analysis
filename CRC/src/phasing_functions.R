extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

lm_eqn <- function(df, x, y){
  m <- lm(y ~ 0 + x, df);
  eq <- substitute(italic(y) == a %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

#' Compute the tumour methylation rate 
#' 
#' \code{compute_tumour_methylome} returns the df with the tumour methylome
#' information added for each CpG
#'
#' @param df_mb_mn dataframe with each CpG and their coverage, counts methylated, 
#' methylation rate and copy number and matched normal methylation info
#' @param p Numerical - Sample purity estimates
#' 
#' @return A dataframe for each sample_id with the tumour methylome added
compute_tumour_methylome <- function (df_mb_mn, p, min_cov_t = 3, gender, build="hg19") {
  flag <- NULL
  
  # convert factors to characters
  p <- as.numeric(p)
  df_mb_mn$CN <- as.numeric(as.character(df_mb_mn$CN))
  
  # Get normal copy number on sex chrom. X in males has CN = 1 outside par regions
  CN_Norm <- rep(2, nrow(df_mb_mn))
  if(gender=="XY"){
    if(build=="hg19"){
      flag <- 
        df_mb_mn[, ifelse(chrom == "X" & start > 60001 & end < 2699520, TRUE, # PAR1
                   ifelse(chrom == "Y" & start > 10001 & end < 2649520, TRUE, # PAR1 
                   ifelse(chrom == "Y" & start > 59034050 & end < 59363566, TRUE, # PAR2
                   ifelse(chrom == "X" & start > 154931044 & end < 155260560, TRUE, # PAR2
                        FALSE))))]
    } 
    if(build=="hg38"){
      flag <- 
        df_mb_mn[, ifelse((chrom == "X" | chrom == "Y") & start > 10001 & end < 2781479, TRUE,
                   ifelse(chrom == "Y" & start > 56887903 & end < 57217415, TRUE,
                   ifelse(chrom == "X" & start > 155701383 & end < 156030895, TRUE, 
                          FALSE)))]
    }
    flag <- flag == FALSE & (df_mb_mn$chrom %in% c("X", "Y"))
    CN_Norm[flag] <- 1
  }
  print(table(CN_Norm))
  
  # get SNV and WT tumour copy number
  df_mb_mn[, "CN_n" := CN_Norm]
  df_mb_mn[, "CN_2" := df_mb_mn[, .(x = seq(1, CN, 1)[which.min(
    abs(seq(1, CN, 1) - SNV_mcn_methylation))]), by = .(1:nrow(df_mb_mn))]$x]
  df_mb_mn[, "CN_1" := CN - CN_2]
  
  # Compute averages per allele  
  df_mb_mn[, m_b := ifelse(is.na(m_1) & is.na(m_2), NA,
                    ifelse(is.na(m_1), M_2 / total_counts_m_2, 
                    ifelse(is.na(m_2), M_1 / total_counts_m_1,  
                          (M_1 + M_2) / (total_counts_m_1 + total_counts_m_2))))]
  
  # Compute allele-specific value  
  df_mb_mn[, m_t_2 := ifelse(is.na(m_2) | is.na(m_n), NA,
                      ifelse(CN_2 == 0, 0, M_2 / total_counts_m_2))]
  df_mb_mn[, m_t_1 := ifelse(is.na(m_1) | is.na(m_n), NA,
                      ifelse(CN_1 == 0, 0,
                            (m_b * (p * CN + (1-p) * CN_Norm) -
                             m_n  * (1-p) * CN_Norm -
                             p * m_t_2 * CN_2) / (p * CN_1)))]
  df_mb_mn[, m_t := ifelse(is.na(m_n), NA, 
                    ifelse(is.na(m_1)|CN_1==0,  
                          (m_b * (p * CN_2 + (1-p) * CN_Norm) - 
                           m_n  * (1-p) * CN_Norm) / (p * CN_2), 
                    ifelse(is.na(m_2)|CN_2==0, 
                          (m_b * (p * CN_1 + (1-p) * CN_Norm) - 
                           m_n  * (1-p) * CN_Norm) / (p * CN_1), 
                          (m_b * (p * (CN) + (1-p) * CN_Norm) - 
                           m_n  * (1-p) * CN_Norm) / (p * (CN)))))]
  
  if(!is.null(flag) & sum(flag == TRUE)>0){
    df_mb_mn$m_t_1[flag] <- NA
    df_mb_mn$m_t_2[flag] <- NA
    df_mb_mn$m_t[flag] <- NA
  }
  
  # Remove NAs (should not be any)
  df_mb_mn <- df_mb_mn[!is.na(m_t_1)|!is.na(m_t_2),]
  
  return(df_mb_mn)
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

HDIofICDF = function(ICDFname, credMass=0.99, tol=1e-4, ... ) {
  
  incredMass = 1.0 - credMass
  
  optInfo = optimize(f = intervalWidth, interval = c(0, incredMass), 
                     ICDFname=ICDFname , credMass=credMass , tol=tol , ... )
  
  HDIlowTailPr = optInfo$minimum
  vec <- setNames(object = ICDFname(c(HDIlowTailPr, credMass+HDIlowTailPr), ... ), 
                  nm = c("low", "high"))
  return(vec)
}