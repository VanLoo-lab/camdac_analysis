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


#' Compute the tumour methylation rate 
#' 
#' \code{get_phased_pure_tumour_methylome} returns the df with the tumour methylome information added for each CpG
#'
#' @param df_mb_mn dataframe with each CpG and their coverage, counts methylated, 
#' methylation rate and copy number and matched normal methylation info
#' @param p Numerical - Sample purity estimates
#' @param sex Character - Either "XX" or "XY" for females and males respectively
#' 
#' @return A dataframe for each sample_id with the tumour methylome added
get_phased_pure_tumour_methylome <- function (df_mb_mn , p, sex) {
  
  # convert factors to characters
  p <- as.numeric(p)
  df_mb_mn$CN <- as.numeric(as.character(df_mb_mn$CN))
  
  # Get allele specific copy number in balanced regions
  df_mb_mn[nA==nB, c("n1", "n2") := list(nA, nA)]
  
  # Get allele specific copy number with phasing info
  df_mb_mn[nA!=nB&phasing!="unknown",n1:=ifelse((phased_seg_BAF<0.5&genotype_1=="ref")|(phased_seg_BAF>0.5&genotype_1=="alt"), nA, nB)]
  df_mb_mn[nA!=nB&phasing!="unknown",n2:=ifelse((phased_seg_BAF<0.5&genotype_1=="ref")|(phased_seg_BAF>0.5&genotype_1=="alt"), nB, nA)]
  
  # Get allele specific copy number w/o phasing info
  df_mb_mn[nA!=nB&phasing=="unknown", n1:=ifelse((POS_genotype_BAF_t<0.5&genotype_1=="ref")|(POS_genotype_BAF_t>0.5&genotype_1=="alt"), nA, nB)]
  df_mb_mn[nA!=nB&phasing=="unknown", n2:=ifelse((POS_genotype_BAF_t<0.5&genotype_1=="ref")|(POS_genotype_BAF_t>0.5&genotype_1=="alt"), nB, nA)]
  
  # Flag CG-altering SNPs
  df_mb_mn[,n1:=ifelse(is.na(CpG_BAF_n), n1, 
                       ifelse((label=="CG-forming"&genotype_1=="alt")|(label=="CG-destroying"&genotype_1=="ref"),
                              n1,0))]

  df_mb_mn[,n2:=ifelse(is.na(CpG_BAF_n), n2, 
                       ifelse((label=="CG-forming"&genotype_2=="alt")|(label=="CG-destroying"&genotype_2=="ref"),
                              n2,0))]
  
  # Get normal copy number (CN = 1 per allele or 0 for the CG-destroying alleles at CG SNPs)
  df_mb_mn[, n1_n := ifelse(is.na(CpG_BAF_n), 1, 
                            ifelse((label=="CG-forming"&genotype_1=="alt")|(label=="CG-destroying"&genotype_1=="ref"),1,0))]
  df_mb_mn[, n2_n := ifelse(is.na(CpG_BAF_n), 1, 
                            ifelse((label=="CG-forming"&genotype_2=="alt")|(label=="CG-destroying"&genotype_2=="ref"),1,0))]
  
  # Get normal copy number on sex chrom. X in males has CN = 1 outside par regions
  if(sex=="XY"){
    if(build=="hg19"){
      flag <- ifelse(df_mb_mn$chrom=="X" & df_mb_mn$start>60001 & df_mb_mn$end<2699520, TRUE, # PAR1
                     ifelse(df_mb_mn$chrom=="Y" & df_mb_mn$start>10001 & df_mb_mn$end<2649520, TRUE, # PAR1     
                            ifelse(df_mb_mn$chrom=="Y" & df_mb_mn$start>59034050 & df_mb_mn$end<59363566, TRUE, # PAR2
                                   ifelse(df_mb_mn$chrom=="X" & df_mb_mn$start>154931044 & df_mb_mn$end<155260560, TRUE, # PAR2
                                          FALSE))))
    } 
    if(build=="hg38"){
      flag <- ifelse((df_mb_mn$chrom=="X" | df_mb_mn$chrom=="Y")& df_mb_mn$start>10001 & df_mb_mn$end<2781479,TRUE,
                     ifelse(df_mb_mn$chrom=="Y" & df_mb_mn$start>56887903 & df_mb_mn$end<57217415, TRUE,
                            ifelse(df_mb_mn$chrom=="X" & df_mb_mn$start>155701383 & df_mb_mn$end<156030895, TRUE, 
                                   FALSE)))
      
    }
    flag <- flag==FALSE & (df_mb_mn$chrom%in%c("X", "Y")) 
    df_mb_mn <- df_mb_mn[!flag]
  }
  
  
  # Compute m_t  
  df_mb_mn[, m_t_1 := ifelse(n1==0, as.numeric(NA), 
                             (m_1 * (p * n1  + (1-p) * n1_n) - m_n_1  * (1-p) * n1_n) / 
                               (p * n1))]
  df_mb_mn[, m_t_2 := ifelse(n2==0, as.numeric(NA),
                             (m_2 * (p * n2  + (1-p) * n2_n) - m_n_2  * (1-p) * n2_n) / 
                               (p * n2))]
  
  # Remove NAs (should not be any)
  df_mb_mn <- df_mb_mn[!is.na(df_mb_mn$m_t_1)|!is.na(df_mb_mn$m_t_2),]
  
  return(df_mb_mn)
  
}

#' Get DMPs
#' 
#' \code{get_DMPs} returns the df with annotated q-value for each CpG
#'
#' @param path Path to bulk methylation RData object exclude filename (/path/to/file/)
#' @param df Bulk methylation object filename with annotated normal info
#' @param n_cores Number of cores to do the statistical testing over
#' 
#' @return A GRanges object with all the CpG loci, their coverage, counts methylated and methylation rate
get_DMPs <- function (alpha_1, alpha_2, beta_1, beta_2,n_cores) { 
  
  # Evan Miller's closed form solution for the probability that a draw from a beta dist is greater than anoter, 
  # p_A ~ Beta(alpha_A, beta_A), p_B ~ Beta(alpha_B, beta_B). In this case, A is the normal and B is the bulk tumour. 
  # Alpha is the counts methylated, beta is counts unmethylated
  # this equation output the probability that p_n > p_b
  h <- function(alpha_n, beta_n, alpha_b, beta_b){
    j <- seq.int(0, round(alpha_b)-1)
    log_vals <- (lbeta(alpha_n + j, beta_n + beta_b) - log(beta_b + j) -
                   lbeta(1 + j, beta_b) - lbeta(alpha_n, beta_n))
    1 - sum(exp(log_vals))
  }
  
  # set up results object
  n<-length(alpha_1)
  result <- cbind(numeric(n))
  
  #Define Variables
  alpha_n=alpha_1+1
  beta_n=beta_1+1
  alpha_b=alpha_2+1
  beta_b=beta_2+1
  # Although a B(0.5, 0.5) prior is more suited, to methylation rates, alpha_b must be an integer value for h() to hold.
  # We therefore add a pseudocount of 1 which is equivalent to B(1, 1) as alpha and beta values go to zero.
  
  idx <- which(rowSums(!is.na(cbind(alpha_n,beta_n,alpha_b,beta_b)))==4)
  m_n <-  alpha_n /(alpha_n + beta_n)
  m_b <-  alpha_b /(alpha_b + beta_b)
  
  result[idx,1] <- 
    mcmapply(
      function(alpha_n, beta_n, alpha_b, beta_b, m_n, m_b) {
        prob <- NULL
        if(m_b <= m_n){
          prob <- h(alpha_n = alpha_n, beta_n=beta_n, alpha_b=alpha_b, beta_b=beta_b)
        } else {
          prob <- h(alpha_n = alpha_b, beta_n=beta_b, alpha_b=alpha_n, beta_b=beta_n)  
        }
        if(is.null(prob)){prob <- as.numeric(NA)}
        return(prob)
      },
      alpha_n=alpha_n[idx],beta_n=beta_n[idx],alpha_b=alpha_b[idx],beta_b=beta_b[idx],
      m_n=m_n[idx],m_b=m_b[idx], mc.cores=n_cores)
    
  result[!idx,1] <- as.numeric(NA)
  prob_DMP <- result[,1]
  prob_DMP <- ifelse(is.na(prob_DMP), as.numeric(NA),
                     ifelse(result[idx,1]>1, 1, 
                            ifelse(result[idx,1]<0, 0, 
                                   result[idx,1])))
  
  DMP_counts_stats_sig <- sum(prob_DMP>0.99, na.rm=TRUE)
  cat(paste0("DMP_counts :",DMP_counts_stats_sig,"\n")) 
  cat(paste0("DMP_prop :", DMP_counts_stats_sig/n,"\n")) 
  # overcalls as expect, Fisher and z-test also overcall
  
  return(prob_DMP)
  
  # clean-up
  rm(result,n,alpha_n,beta_n,alpha_b,beta_b, m_n, m_b)
} 

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
#' is a scalar between 0 and 1, indicating the mass within the credible
#' interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
HDIofMCMC = function(M_b,UM_b,M_n,UM_n,p,CN,CN_n,credMass=0.99) {
  
  # create blank object
  HDIlim <- c(lower = as.numeric(NA), 
              upper = as.numeric(NA))
  
  # sampleVec is a vector of representative values for the methylation rate 
  # probability distribution at a CpG
  if(sum(is.na(c(M_b, UM_b, M_n, UM_n, CN, CN_n)))==0 & CN>0){
    sampleVec1 <- rbeta(n = 2000, shape1 = M_b+1, shape2 = UM_b+1)
    sampleVec2 <- rbeta(n = 2000, shape1 = M_n+1, shape2 = UM_n+1)
    
    # get the mb and mn constant
    cons1 <- p*CN + CN_n*(1-p)
    cons2 <- CN_n*(1-p)
    
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
  } 
  
  return(HDIlim)
}

# Calculate HDI overlaps by simulation
#
# Computes highest density interval from a sample of representative values,
#   estimated as shortest credible interval for a unimodal distribution
# Arguments:
#' @param M_b_1 vector with counts methylated for the tumour allele 1
#' @param UM_b_1 vector with counts unmethylated for the tumour allele 1
#' @param M_n_1 vector with counts methylated for the tumour allele 1
#' @param UM_n_1 vector with counts unmethylated for the tumour allele 1
#' @param CN_1 tumour CpG allele 1 copy number 
#' @param CN_n_1 normal CpG allele 1 copy number 
#' @param M_b_2 vector with counts methylated in the tumour sample (or allele) 2
#' @param UM_b_2 vector with counts unmethylated in the tumour sample (or allele) 2
#' @param M_n_2 vector with counts methylated in the normal sample (or allele) 2
#' @param UM_n_2 vector with counts unmethylated in the normal sample (or allele) 2
#' @param CN_2 tumour CpG allele 2 copy number 
#' @param CN_n_2 normal CpG allele 2 copy number 
#' @param p tumour purity
#' is a scalar between 0 and 1, indicating the mass within the credible
#' interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
MCMC_overlap_pure = function(M_b_1,UM_b_1,M_n_1,UM_n_1,CN_1,CN_n_1,
                             M_b_2,UM_b_2,M_n_2,UM_n_2,CN_2,CN_n_2, p) {
  # M_b_1 <- 80; UM_b_1 <- 20; CN_1 <- 2 
  # M_n_1 <- 45; UM_n_1 <- 55; CN_n_1 <- 1
  # M_b_2 <- 20; UM_b_2 <- 80; CN_2 <- 1
  # M_n_2 <- 45; UM_n_2 <- 55; CN_n_2 <- 1
  
  # create blank object
  ov <- as.numeric(NA)
  
  # sampleVec is a vector of representative values for the methylation rate 
  # probability distribution at a CpG
  if(sum(is.na(c(M_b_1, UM_b_1, M_n_1, UM_n_1, CN_1, CN_n_1,
                 M_b_2, UM_b_2, M_n_2, UM_n_2, CN_2, CN_n_2)))==0 & 
     CN_1>0 & CN_2>0 & CN_n_1>0 & CN_n_2>0){
    for(suffix in c(1, 2)){
       M_b <- unlist(mget(paste0("M_b_", suffix)))
       UM_b <- unlist(mget(paste0("UM_b_", suffix)))
       CN <- unlist(mget(paste0("CN_", suffix)))
       M_n <- unlist(mget(paste0("M_n_", suffix)))
       UM_n <- unlist(mget(paste0("UM_n_", suffix)))
       CN_n <- unlist(mget(paste0("CN_n_", suffix)))
     
        sampleVec1 <- rbeta(n = 2000, shape1 = M_b+1, shape2 = UM_b+1)
        sampleVec2 <- rbeta(n = 2000, shape1 = M_n+1, shape2 = UM_n+1)
        
        # get the mb and mn constant
        cons1 <- p*CN + CN_n*(1-p)
        cons2 <- CN_n*(1-p)
        
        # obtain the mt sample vector
        sampleVec_scaled <-  (cons1*sampleVec1 - cons2*sampleVec2)/(p*CN)
        
        # sort the mt sample vector
        sortedPts <- sort(sampleVec_scaled)
        assign(paste0("sortedPts_", suffix), sortedPts)
        
        # get the mean methylation rate
        m <- mean(sortedPts, na.rm=TRUE)
        assign(paste0("m_", suffix), m)
      } 
  
   if(min(sortedPts_1)>max(sortedPts_2)|min(sortedPts_2)>max(sortedPts_1)){
       ov <- 0
   } else {
       # boundaries <- list(from = c(0, 0), to = c(1,1))
       # # library(overlapping)
       # ov <- overlap(list(sortedPts_1, sortedPts_2), nbins = 2000, plot = FALSE, partial.plot = FALSE,  boundaries= boundaries)$OV
       if(m_1>m_2){
         n <- which.max(sortedPts_1<sortedPts_2)
         ov <- n/2000
       }
       if(m_1<m_2){
         n <- which.max(sortedPts_2<sortedPts_1)
         ov <- n/2000
       }
    }
  }
  return(ov)
}

# Calculate HDI overlaps by simulation
#
# Computes highest density interval from a sample of representative values,
#   estimated as shortest credible interval for a unimodal distribution
# Arguments:
#' @param M_1 vector with counts methylated for the tumour allele 1
#' @param UM_1 vector with counts unmethylated for the tumour allele 1
#' @param CN_1 tumour CpG allele 1 copy number 
#' @param M_2 vector with counts methylated in the tumour sample (or allele) 2
#' @param UM_2 vector with counts unmethylated in the tumour sample (or allele) 2
#' @param CN_2 tumour CpG allele 2 copy number 
#' is a scalar between 0 and 1, indicating the mass within the credible
#' interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
MCMC_overlap_beta = function(M_1,UM_1,M_2,UM_2,BAF) {
  # M_1 <- 40; UM_1 <- 10; CN_1 <- 1
  # M_2 <- 20; UM_2 <- 30; CN_2 <- 1

  # create blank object
  ov <- as.numeric(NA)
  
  # sampleVec is a vector of representative values for the methylation rate 
  # probability distribution at a CpG
  if(sum(is.na(c(M_1,UM_1,M_2,UM_2)))==0){  
      if(sum(c(M_1,UM_1))>0 & sum(c(M_2,UM_2))>0){ 
        if(is.na(BAF) | (!is.na(BAF) & (BAF<0.3 | BAF>0.7))){
            # get vectors
            sampleVec1 <- rbeta(n = 2000, shape1 = M_1+1, shape2 = UM_1+1)
            sampleVec2 <- rbeta(n = 2000, shape1 = M_2+1, shape2 = UM_2+1)
            
            # sort the mt sample vector
            sortedPts_1 <- sort(sampleVec1)
            m_1 <- mean(sortedPts_1, na.rm=TRUE)
            sortedPts_2 <- sort(sampleVec2)
            m_2 <- mean(sortedPts_2, na.rm=TRUE)
            
            if(min(sortedPts_1)>max(sortedPts_2)|min(sortedPts_2)>max(sortedPts_1)){
              ov <- 0
            } else {
              # boundaries <- list(from = c(0, 0), to = c(1,1))
              # ov <- overlap(list(sortedPts_1, sortedPts_2), nbins = 2000, plot = FALSE, 
              #               partial.plot = FALSE,  boundaries= boundaries)$OV
              # ov <- unname(ov)
              # }
              if(m_1>m_2){
                n <- which.max(sortedPts_1<sortedPts_2)
                ov <- n/2000
              }
              if(m_1<m_2){
                n <- which.max(sortedPts_2<sortedPts_1)
                ov <- n/2000
              }
            }
        }
      }
  }
  return(ov)
}


# Arguments:
#'  @param i 
#'  @param j 
#'  @param n_cores
#'  @param bam_file
#'  @param loci_subset
#'  @param df_POS
#'  @param test
#'  @param mq

# function to get SNP and CpG allele counts
get_reads <- function(i, j = NULL, n_cores = NULL, bam_file, loci_subset, df_POS, test, mq=0){
  
  # Store segements file in tmp object
  loci_subset_2 <- loci_subset
  
  # If this is a test run - downsample
  if(test==TRUE){
    set.seed(0)
    loci_subset_2 <- sort(loci_subset_2[sample(1:length(loci_subset_2), 5E3)])
  }
  
  # If running in parallel, further subset the genomic segments file (input for scanBam)
  if(!is.null(j)){
    loci_subset_2$coreID <- as.numeric(ceiling(1:length(start(loci_subset_2))/(length(start(loci_subset_2))/n_cores+1)))
    loci_subset_2 <- loci_subset_2[loci_subset_2$coreID==j,]
  }
  
  # remove metadata
  mcols(loci_subset_2) <- NULL
  
  # Set what argument to speed up scanBam
  whats <- c("qname", "rname", "strand", "pos", "qwidth", "mapq", "seq", "qual") # scanBamWhat()
  idxFile <- paste(bam_file, ".bai", sep = "")
  
  # Set info that we want to obtain for each position
  bamFile <- BamFile(bam_file)
  params <- Rsamtools::ScanBamParam(which=loci_subset_2, what=whats)
  
  # Run ScanBam 
  #bam <- Rsamtools::scanBam(file=bam_file, index=idxFile, param=params)
  bam <- Rsamtools::scanBam(file=bamFile, param=params)
  rm(params, loci_subset_2)
  
  # Remove entries with zero read
  bam_filtered <- bam[which(sapply(bam, function(x) length(x[["qname"]])) != 0)]
  rm(bam)
  
  # Unlist bam
  df <- data.table(dplyr::bind_rows(lapply(bam_filtered, as.data.frame.list, stringsAsFactors=FALSE)))
  # duplicates created when unlisting and due to overlapping reads, to be removed in later lines
  rm(bam_filtered)
  
  # Format pileup df for transformation into GRanges Object
  df$end <- df$pos + (df$qwidth - 1)
  colnames(df)[c(which(colnames(df) == "rname"),which(colnames(df) == "pos"))] <- c("chrom", "start")
  rownames(df) <- c(1:nrow(df))
  
  cols <- c("qname", "chrom", "strand", "start", "end", "mapq", "seq", "qual")
  Sys.sleep(1)
  
  df_bam <- df[,cols,with=FALSE] 
  rm(df)
  
  # The pileup unlisting produces duplicate rows where more than one loci overlaps with the same molecule
  x <- do.call(paste, c(df_bam[,cols,with=FALSE], sep="-")) 
  df_bam_dedup <- df_bam[!duplicated(x),]
  rm(df_bam, cols, x)
  
  # Create Granges - Each row represents a single read 
  gr_bam <- makeGRangesFromDataFrame(df_bam_dedup, keep.extra.columns = TRUE)
  rm(df_bam_dedup)
  
  # Change seqnames if necessary
  if(UCSC == FALSE){seqlevelsStyle(loci_subset) <- "Ensembl"}
  
  # Extract info per CpG/SNP loci
  overlaps <- mergeByOverlaps(gr_bam, loci_subset)
  rm(gr_bam)
  
  # store as data.table
  df_pileup <- data.table(qname=as.character(overlaps$qname), strand=as.character(strand(overlaps[,"gr_bam"])),
                          chrom=as.character(seqnames(overlaps[,"gr_bam"])),read.start=start(overlaps[,"gr_bam"]), 
                          read.end=end(overlaps[,"gr_bam"]),POS=overlaps$POS, width = width(overlaps[,"loci_subset"]),
                          start=start(overlaps[,"loci_subset"]),end=end(overlaps[,"loci_subset"]),ref=overlaps$ref,
                          alt=overlaps$alt,seq=overlaps$seq,qual=overlaps$qual,mq=overlaps$mapq,stringsAsFactors = FALSE)
  #rm(overlaps)
  
  cols<-c("qname","chrom","POS","width","start","ref","alt")
  x <- do.call(paste, c(df_pileup[,cols,with=FALSE], sep="-")) 
  df_pileup_dedup <- df_pileup[!duplicated(x),]
  rm(df_pileup, cols, x)
  
  # Flag 5'CCGG and rm 3' CCGG but keep any SNPs. 
  df_pileup_dedup[, CCGG := (ifelse(width != 4 , "non-CCGG",
                                    ifelse(((strand == "+" & read.start < end & read.start >= start) | 
                                              (strand == "-" & read.end >= start & read.end <= end)), "5pCCGG",
                                           ifelse(((strand == "-" & read.start <= end & read.start >= start) |
                                                     (strand == "+" & read.end >= start & read.end <= end)), "3pCCGG",
                                                  ifelse(!is.na(POS), "CCGG-SNP", "CCGG-SNV")))))]
  
  # Assuming complete digestion, non-3p/5p CCGG can occur due to: 
  # a germline/tumour mutation destroying CCGG motif (CCGG-SNV) or
  # a SNP allele destroying the CCGG motif (CCGG-SNP) 
  
  # Remove mid-read CCGGs due to incomplete digestion or reference-mismtach (germline or somatic variant but keep SNPs)
  df_pileup_dedup[, flag := ifelse(CCGG=="CCGG-SNV",FALSE,TRUE)]
  df_pileup_dedup <- df_pileup_dedup[flag == TRUE,] 
  df_pileup_dedup[, flag := NULL] 
  
  # Flag SNPs at CCGG-destroying SNPs 
  df_pileup_dedup[, flag := (width==4 & !is.na(POS))] 
  df_pileup_dedup[, CCGG_SNP := (ifelse(flag == FALSE, "n/a",
                                        ifelse(CCGG=="5pCCGG" & # Flag 5' SNPs inside of read 
                                                 ((strand=="+" & POS!=start) | (strand=="-" & POS!=end)), "yes",
                                               ifelse(CCGG=="3pCCGG" & # Flag 3' SNPs inside of read 
                                                        ((strand=="+" & POS==start) | (strand=="-" & POS==end)), "yes", 
                                                      ifelse(CCGG=="CCGG-SNP","yes","no")))))]
  df_pileup_dedup[, flag := NULL]
  
  # Extract 5'CCGG alleles
  df_pileup_dedup[, flag := (CCGG=="5pCCGG")]
  df_pileup_dedup[, tmp := (start-read.start)]
  df_pileup_dedup[, alleles.CCGG := ifelse(flag == FALSE, NA, paste(substr(as.character(seq),tmp+1,tmp+4),strand,sep=""))]
  df_pileup_dedup[, flag := NULL]
  
  # Annotate allowed 5' alleles CGG+, TGG+, CCG- or CCA- 
  z <- grep("(^CC([AG]{1})([ACGNT]{0,1})([-])$)|(([ACGNT]{0,1})([CT]{1})GG([+]))",df_pileup_dedup$alleles.CCGG)
  
  # Remove reads where 5' CCGG isn'TRUE match [CT]CGG+, [CT]TGG+, CCG[GA]- or CCA[GA]- 
  # These are usually due to poor alignment (i.e. tend to be C rich )
  df_pileup_dedup[, flag := logical(length=nrow(df_pileup_dedup))]
  df_pileup_dedup[z, flag := TRUE] ; rm(z)
  
  # Find all misaligned reads
  qnames <- unique(df_pileup_dedup[,(qname[flag==FALSE & CCGG=="5pCCGG"])])
  
  # Keep reads which have an alternative overlapping 5'CCGG w/ suitable alleles 
  tmp2 <- unique(df_pileup_dedup[,(qname[qname%in%qnames & flag==TRUE])])
  qnames <- qnames[!qnames%in%tmp2];rm(tmp2)
  
  # Remove misaligned reads
  df_pileup_dedup<-df_pileup_dedup[!qname%in%qnames,];rm(qnames)
  
  # Rename remaining non-5p CCGGs which are due to SNP / SNV
  flag2 <- df_pileup_dedup[, (flag==FALSE & CCGG=="5pCCGG")]
  df_pileup_dedup[flag2 & is.na(POS), CCGG := "CCGG-SNV"] 
  df_pileup_dedup[flag2 & !is.na(POS), CCGG := "CCGG-SNP"] 
  # Note. SNPs/SNVs are unlikely to form a valid 5p CCGG.
  # To do so, it would have to turn a motif overlappingto the ref CCGG to another CCGG 
  # (i.e. -CCGGCGG- to -CCGCCGG-)
  
  # clean-up
  df_pileup_dedup[, flag := NULL] ; rm(flag2)
  
  # # print diganostics
  # cat("Misaligned reads removed", "Adapter contamination removed!", sep = "\n", file = file_stdout_2, append = TRUE)
  
  # Extract CpG and SNP alleles
  df_pileup_dedup[, offset := ifelse(width == 4, 2, 1)]
  df_pileup_dedup[, alleles.dinucs :=  ifelse(width>=2&CCGG!="3pCCGG", 
                                              paste(substr(seq, tmp+offset, tmp+offset+1),strand, sep=""), NA)]
  df_pileup_dedup[, tmp2 := (POS-read.start+1)]
  df_pileup_dedup[, alleles.SNP := ifelse(!is.na(POS)&CCGG_SNP!="no", paste(substr(seq,tmp2,tmp2),strand, sep=""), NA)]
  
  # Extract CpG and SNP quality from the base quality string
  df_pileup_dedup[, qual.dinucs := ifelse(width>=2&CCGG!="3pCCGG",substr(qual,tmp+offset, tmp+offset+1), NA)]
  df_pileup_dedup[, qual.SNP := ifelse(!is.na(POS)&CCGG_SNP!="no",substr(qual,tmp2,tmp2),NA)]
  
  # get read length
  df_pileup_dedup[, read_length := (read.end-read.start+1)]
  
  # Concatenate results
  cols <- c("chrom", "strand", "start", "end", "width", "CCGG", "alleles.CCGG", "POS", "ref", "alt",
            "alleles.dinucs", "alleles.SNP", "qual.dinucs", "qual.SNP", "mq", "qname", "read_length")
  assign("df_pileup_alleles",df_pileup_dedup[,cols,with=FALSE]);rm(df_pileup_dedup,cols)
  
  # Flag SNPs and CpGs where BQ < 20 
  df_pileup_alleles[, flag1 := width >= 2 & grepl(pattern = "([5-9A-K:-@])([5-9A-K:-@])", qual.dinucs)]
  df_pileup_alleles[, flag2 := !is.na(POS) & grepl(pattern = "([5-9A-K:-@])", qual.SNP)]
  df_pileup_alleles[, flag := (flag1==TRUE | flag2==TRUE)] # CpG or SNP has BQ > 20
  
  # Remove positions where neither CG and/or SNP passes the BQ tresholds
  if(sum(df_pileup_alleles$flag)!=nrow(df_pileup_alleles)){
    df_pileup_alleles <- df_pileup_alleles[df_pileup_alleles$flag == TRUE,] }
  df_pileup_alleles[, flag := NULL]
  
  # Remove dinucs at CpGs where only SNPs passes qual treshold, but keep SNP loci for BAF and logR
  df_pileup_alleles[, alleles.dinucs := ifelse(flag1==TRUE, alleles.dinucs, NA)]
  df_pileup_alleles[, qual.dinucs := ifelse(flag2==TRUE, qual.dinucs, NA)]
  df_pileup_alleles[, c("flag1", "flag2") := NULL]
  
  # remove reads where MQ < mq treshold (recommeded is 0)
  if(mq > 0){df_pileup_alleles <- df_pileup_alleles[df_pileup_alleles$mq >= mq,]}
  
  # Assign ref or alt count for each read
  df_pileup_alleles[, Af := ifelse(alleles.SNP == "A+",1,0)]
  df_pileup_alleles[, Ar := ifelse(alleles.SNP == "A-",1,0)]
  df_pileup_alleles[, Cf := ifelse(alleles.SNP == "C+",1,0)]
  df_pileup_alleles[, Cr := ifelse(alleles.SNP == "C-",1,0)]
  df_pileup_alleles[, Gf := ifelse(alleles.SNP == "G+",1,0)]
  df_pileup_alleles[, Gr := ifelse(alleles.SNP == "G-",1,0)]
  df_pileup_alleles[, Tf := ifelse(alleles.SNP == "T+",1,0)]
  df_pileup_alleles[, Tr := ifelse(alleles.SNP == "T-",1,0)]
  
  # Assign ref or alt count for each read
  df_pileup_alleles[, CGr := ifelse(alleles.dinucs == "CG-",1,0)]
  df_pileup_alleles[, CGf := ifelse(alleles.dinucs == "CG+",1,0)]
  df_pileup_alleles[, CAr := ifelse(alleles.dinucs == "CA-",1,0)]
  df_pileup_alleles[, TGf := ifelse(alleles.dinucs == "TG+",1,0)]
  
  # Assign fragments breakpoint
  df_pileup_alleles[, CCGG := ifelse(CCGG == "5pCCGG",1,0)]
  
  # Create tmp object for alt and ref calculations
  tmp <- df_pileup_alleles[!is.na(POS),.(chrom, strand, width, POS, ref, alt, mq, qname, Af, Ar, Cf, Cr, Gf, Gr, Tf, Tr)]
  
  # only take het SNPs (hom SNP neighbouring het SNPs will have been piled - up)
  df_POS$flag <- "Het"
  if(UCSC==T&df_POS$chrom[1]!="chr1"){
    df_POS[, chrom2 := chrom]
    df_POS[, chrom := paste0("chr", chrom2)]
    df_POS[, chrom2 := NULL]
  }
  tmp <- merge(tmp, df_POS, by = c("chrom", "POS", "ref", "alt"))
  
  # Get SNP in ref:alt format 
  tmp[, SNP := (ifelse(!is.na(ref), paste(ref, alt, sep = ""), NA))]
  tmp[, offset := (ifelse(width == 4, 1, 0))]
  
  # Calculate the ref counts
  tmp[,ref_counts := (ifelse(is.na(ref), NA,
                      ifelse(grepl(pattern="AC",SNP), Af + Ar,
                      ifelse(grepl(pattern="CA",SNP), Tf + Cr + Cf,
                      ifelse(grepl(pattern="AG",SNP), Af,
                      ifelse(grepl(pattern="GA",SNP), Gf,      
                      ifelse(grepl(pattern="AT",SNP), Af + Ar,
                      ifelse(grepl(pattern="TA",SNP), Tf + Tr,
                      ifelse(grepl(pattern="GT",SNP), Gf + Ar + Gr,
                      ifelse(grepl(pattern="TG",SNP), Tf + Tr,
                      ifelse(grepl(pattern="CG",SNP), Tf + Cr + Cf,
                      ifelse(grepl(pattern="GC",SNP), Gf + Ar + Gr,
                      ifelse(grepl(pattern="CT",SNP), Cr,
                      ifelse(grepl(pattern="TC",SNP), Tr, NA))))))))))))))]
  # We consider all bisulfite converted alleles and that C methylation is possible both in and outside CpG context
  
  # Calculate the alt counts
  tmp[, alt_counts := (ifelse(is.na(ref), NA,
                       ifelse(grepl(pattern="AC",SNP), Tf + Cr + Cf,
                       ifelse(grepl(pattern="CA",SNP), Af + Ar,
                       ifelse(grepl(pattern="AG",SNP), Gf,
                       ifelse(grepl(pattern="GA",SNP), Af,      
                       ifelse(grepl(pattern="AT",SNP), Tf + Tr,
                       ifelse(grepl(pattern="TA",SNP), Af + Ar,
                       ifelse(grepl(pattern="GT",SNP), Tf + Tr,
                       ifelse(grepl(pattern="TG",SNP), Gf + Ar + Gr,
                       ifelse(grepl(pattern="CG",SNP), Gf + Ar + Gr,
                       ifelse(grepl(pattern="GC",SNP), Tf + Cr + Cf,
                       ifelse(grepl(pattern="CT",SNP), Tr,
                       ifelse(grepl(pattern="TC",SNP), Cr, NA))))))))))))))]  
  
  # get the genotype for each row (each SNP position)
  tmp <- tmp[,.("ref_counts" = sum(ref_counts, na.rm = TRUE), "alt_counts" = sum(alt_counts, na.rm = TRUE)), 
             keyby = .(qname,strand,chrom,POS,ref,alt)]
  #table(tmp[,(paste(ref_counts,alt_counts))]) 
  # should only have 00,10,01 as each read is either ref (10) or alt (01) or unassigned (00). 
  # duplicated SNP entries lead to ref "20" and alt "02"
  
  # Adjust genotypes at SNP which have duplicated row entries 
  # For example, this can occur due to CCG>CGG SNP (one enrty for each potential CpG)
  tmp$alt_counts <- tmp[,(ifelse(alt_counts==2, 1, alt_counts))]
  tmp$ref_counts <- tmp[,(ifelse(ref_counts==2, 1, ref_counts))]
  tmp$genotype <- tmp[,(ifelse(alt_counts==1&ref_counts==0, "alt", 
                               ifelse(ref_counts==1&alt_counts==0, "ref", NA)))]
  
  # remove reads with unknow genotype
  tmp <- tmp[!is.na(tmp$genotype),]
  
  # store SNP ID
  tmp[,SNP_ID := paste(POS,ref,alt,sep=":")]
  
  # simplify the genotype for each read
  tmp2 <- tmp[order(POS),.(num_mismatches = length(unique(SNP_ID[genotype=="alt"])),
                           genotype = ifelse(sum(grepl("ref",genotype))==length(genotype),"ref", 
                                             ifelse(sum(grepl("alt",genotype))==length(genotype),"alt", 
                                                    paste0(genotype,collapse = ","))),
                           POS_genotype = paste0(sort(POS),collapse = ","), 
                           num_POS_genotype = length(unique(SNP_ID))),
              by=.(qname,strand)]
  
  
  # set matching ids 
  tmp2[, ID := paste(qname,strand,sep=";")]
  df_pileup_alleles[, ID := paste(qname,strand,sep=";")]
  
  # Annotate genotype and keep track of reads with > 1 SNP 
  df_pileup_alleles[, genotype := tmp2$genotype[match(ID, tmp2$ID)]]
  df_pileup_alleles[, POS_genotype := tmp2$POS_genotype[match(ID, tmp2$ID)]]
  df_pileup_alleles[, num_POS_genotype := tmp2$num_POS_genotype[match(ID, tmp2$ID)]]
  df_pileup_alleles[, num_mismatches := tmp2$num_mismatches[match(ID, tmp2$ID)]]
  #df_pileup_alleles[, mix_haplotype := tmp2$mix_haplotype[match(ID, tmp2$ID)]]
  rm(tmp2)
  
  # Remove positions with genotype = NA (i.e. no het SNP on read)
  df_pileup_alleles <- df_pileup_alleles[!is.na(genotype),]
  rm(tmp)
  
  # Flag positions with clear ref/alt genotype
  # ignore V1 >= 3 for now
  flag <- df_pileup_alleles[, .(V1 = length(unique(genotype)), 
                                V2 = list(unique(genotype)), 
                                V3 = list(table(genotype)),
                                V4 = list(unique(POS_genotype))),
                                #V5 = sum(mix_haplotype)),
                            by = .(chrom, start, end, width, POS, ref, alt)]
  flag <- flag[V1==2,]
  
  # Keep position with clear ref / alt alleles
  df_pileup_alleles <- merge(df_pileup_alleles, flag[,.(chrom, start, end, width, POS, ref, alt)], 
                             by = c("chrom", "POS", "start", "end", "width","ref", "alt"))
  
  # Calculate total ref or alt count for each position
  df_summary <- df_pileup_alleles[,.("Af" = sum(Af, na.rm = TRUE), "Ar" = sum(Ar, na.rm = TRUE),
                                     "Cf" = sum(Cf, na.rm = TRUE), "Cr" = sum(Cr, na.rm = TRUE),
                                     "Gf" = sum(Gf, na.rm = TRUE), "Gr" = sum(Gr, na.rm = TRUE),
                                     "Tf" = sum(Tf, na.rm = TRUE), "Tr" = sum(Tr, na.rm = TRUE),
                                     "CGf" = sum(CGf, na.rm = TRUE), "CGr" = sum(CGr, na.rm = TRUE),
                                     "TGf" = sum(TGf, na.rm = TRUE), "CAr" = sum(CAr, na.rm = TRUE),
                                     "total_depth" = length(Af), "CCGG" = sum(CCGG, na.rm = TRUE)),#,"mq"=paste0(mq,collapse=",")),
                                  keyby = .(chrom, start, end, width, POS, ref, alt, genotype, POS_genotype)]
  
  # clean-up
  rm(df_pileup_alleles)
  
  #Get SNP in ref:alt format 
  df_summary[, SNP := ifelse(!is.na(ref), paste(ref, alt, sep = ""), NA)]
  df_summary[, offset := ifelse(width == 4, 1, 0)]
  
  #Calculate the ref counts
  df_summary[, ref_counts := ifelse(is.na(ref), NA,
                             ifelse(grepl(pattern="AC",SNP), Af + Ar,
                             ifelse(grepl(pattern="CA",SNP), Tf + Cr + Cf,
                             ifelse(grepl(pattern="AG",SNP), Af,
                             ifelse(grepl(pattern="GA",SNP), Gf,      
                             ifelse(grepl(pattern="AT",SNP), Af + Ar,
                             ifelse(grepl(pattern="TA",SNP), Tf + Tr,
                             ifelse(grepl(pattern="GT",SNP), Gf + Ar + Gr,
                             ifelse(grepl(pattern="TG",SNP), Tf + Tr,
                             ifelse(grepl(pattern="CG",SNP), Tf + Cr + Cf,
                             ifelse(grepl(pattern="GC",SNP), Gf + Ar + Gr,
                             ifelse(grepl(pattern="CT",SNP), Cr,
                             ifelse(grepl(pattern="TC",SNP), Tr, NA)))))))))))))]
  
  #We consider all bisulfite converted alleles and that C methylation is possible both in and outside CpG context
  
  #Calculate the alt counts
  df_summary[, alt_counts := ifelse(is.na(ref), NA,
                             ifelse(grepl(pattern="AC",SNP), Tf + Cr + Cf,
                             ifelse(grepl(pattern="CA",SNP), Af + Ar,
                             ifelse(grepl(pattern="AG",SNP), Gf,
                             ifelse(grepl(pattern="GA",SNP), Af,      
                             ifelse(grepl(pattern="AT",SNP), Tf + Tr,
                             ifelse(grepl(pattern="TA",SNP), Af + Ar,
                             ifelse(grepl(pattern="GT",SNP), Tf + Tr,
                             ifelse(grepl(pattern="TG",SNP), Gf + Ar + Gr,
                             ifelse(grepl(pattern="CG",SNP), Gf + Ar + Gr,
                             ifelse(grepl(pattern="GC",SNP), Tf + Cr + Cf,
                             ifelse(grepl(pattern="CT",SNP), Tr,
                             ifelse(grepl(pattern="TC",SNP), Cr, NA)))))))))))))]          
  
  df_summary[, total_counts := ifelse(is.na(POS), NA, ref_counts+alt_counts)]
  
  # Counts methylated for allele-specific methylation i.e. excluding CG destroying SNPs.
  df_summary[, M := ifelse(width == 1, NA,
                           ifelse(is.na(SNP), CGf + CGr,
                                  ifelse(start+offset == POS & grepl(pattern = "([CT][TC])",SNP), CGr,
                                         ifelse(end-offset == POS &grepl(pattern = "([GA][AG])",SNP), CGf,
                                                CGf + CGr))))]
  
  # Get the total methylation counts --> CG-destroying SNP allele does not contribute. 
  # Will be NA if the CG-destroying allele at CG/CCGG is the only one present
  df_summary[, total_counts_m := ifelse(width == 1, NA,
                                 ifelse(is.na(SNP), TGf + CAr + M,
                                 ifelse(start+offset == POS & grepl(pattern = "([CT][TC])", SNP), CAr + M,
                                 ifelse(end-offset == POS & grepl(pattern = "([GA][AG])", SNP), TGf + M, 
                                        TGf + CAr + M))))]
  
  # Get the counts unmethylated
  df_summary[, UM := total_counts_m-M]
  df_summary[, offset := NULL]
  
  # Calculate the all counts including ones indistinguishable between alt and ref and /or meth and unmethyl
  df_summary[, all_counts := ifelse(is.na(ref), TGf+CAr+CGf+CGr,
                                    ifelse(grepl(pattern="([GA][AG])",SNP),Af+Gf+Ar+Gr,
                                           ifelse(grepl(pattern="([CT][TC])",SNP),Cr+Tr+Tf+Cf,
                                                  total_counts)))]
  
  #Remove rows with no allowed counts
  flag <- df_summary[,is.na(total_counts_m)&is.na(total_counts)]
  if(sum(flag)>0){df_summary <- df_summary[flag==FALSE,]};rm(flag)
  
  # Remove positions with unexpected bases at SNPs (i.e. SNV)
  df_summary[, other_counts :=  total_depth - all_counts]
  flag <- df_summary[, (ifelse(other_counts <= 0.05*total_depth, TRUE, FALSE))]
  if(!is.na(table(flag)["FALSE"])) {
    warning(paste0("Unexpected nucleotides represent > 5% of the alleles at ", sum(table(flag)["FALSE"]), 
                   " positions ( ", round(sum(table(flag)["FALSE"])*100/length(flag), digits = 2),
                   " %) - positions removed.", sep = ""))
    df_summary<- df_summary[flag,]
  };rm(flag)
  
  df <- df_summary[,c("chrom","start","end","width","POS","ref","alt","alt_counts","ref_counts",
                      "total_counts","total_depth","other_counts","all_counts","M","UM","total_counts_m",
                      "Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr","CCGG",
                      "genotype", "POS_genotype"), with = F]
  # clean up
  rm(df_summary)
  colnames(df)[c(1)] <- c("CHR")
  
  # Remove loci with 0 distinguishable ref/alt alleles & no methylation
  df[, flag := ifelse(!is.na(ref)&total_counts>0, TRUE, ifelse(!is.na(total_counts_m)&total_counts_m>=3, TRUE,FALSE))]
  if(!is.na(table(df$flag)["FALSE"])) {
    warning(paste0("Removing loci with zero SNP allele counts and/or insufficient methylation - ", sum(table(df$flag)["FALSE"]), 
                   " positions removed.", sep = ""))
    df <- df[flag == TRUE,]
  }
  df[, flag := NULL]
  
  # Get Ensembl chrom annotation for compatibility with ASCAT
  if(UCSC == TRUE){df[,chrom := as.character(substr(df$CHR, 4, 5))]}
  if(UCSC == FALSE){df[,chrom := CHR]}
  df[,chrom := factor(chrom,levels=c(1:22,"X","Y"),ordered = TRUE)]
  df[,CHR := factor(CHR,levels=paste0("chr",c(1:22,"X","Y")),ordered = TRUE)]
  
  # Compute methylation rate for CpGs
  df[, m := ifelse(width == 1, NA, ifelse(total_counts_m == 0, NA, M / (M + UM)))]
  
  # Compute BAF for SNPs 
  df[, BAF := ifelse(is.na(ref), NA, ifelse(total_counts == 0, NA, alt_counts / total_counts))]
  
  # Format row names and col classes
  rownames(df) <- 1:nrow(df)
  if(getOption("scipen")==0){options(scipen = 999)}
  
  # Ensure there are no NAs before saving
  flag <- df[,!is.na(CHR)]
  if(sum(flag==FALSE)>0){df <- df[flag,]};rm(flag)
  
  # Select columns for saving and put in a logical order
  df <- df[,.SD, .SDcols=c("CHR","chrom","start","end","width","POS","ref","alt",
                           "alt_counts","ref_counts","total_counts","BAF","total_depth","other_counts",
                           "all_counts","M","UM","total_counts_m","m",
                           "Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr","CCGG",
                           "genotype", "POS_genotype")]
  
  # keep rows with both alleles remaining after QC
  flag <- df[,(num_genotypes=length(unique(na.omit(genotype)))), 
             keyby=c("CHR","chrom","start","end","width","POS","ref","alt")]
  flag <- flag[V1==2,]
  df <- merge(df,flag,by=c("CHR","chrom","start","end","width","POS","ref","alt"))
  
  df <- data.frame(df, stringsAsFactors = FALSE);rm(flag)
  return(df)
}