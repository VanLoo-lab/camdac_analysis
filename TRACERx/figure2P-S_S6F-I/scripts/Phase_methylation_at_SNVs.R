#' @description Get methylation phased with SNV
#' @author Lili Larose Cadieux

source("~/.Rprofile")

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
compute_tumour_methylome <-
  function (df_mb_mn, p, min_cov_t = 3, gender, build="hg19") {
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

# Load library
x <- c("GenomicRanges", "Rsamtools", "stringr", "readr", "matrixStats",
       "ggplot2", "gridExtra", "data.table", "parallel", "ggpubr", "fst")
invisible(lapply(x, function(y) 
  suppressWarnings(suppressMessages(
    library(y, character.only = TRUE, quietly = TRUE, 
            warn.conflicts = FALSE)))))
rm(x)

LOH_analysis = FALSE
neo = FALSE

# Set static vars
n_cores=1
build = "hg19"
gender_list <- read_csv("~/rdeocd/txMeth/tx100_RRBS_patients.csv")
gender_list$Gender <- ifelse(gender_list$Gender=="Female", "XX", "XY")
mq <- 0
path_to_CAMDAC <- "~/"
parent_dir <- "~/rdeocd/txMeth/"
bam_dir <- "~/storage/CI_FeberTanic_TracerX/TRACERx/shared/NuGEN_sequencing_data/"

# load reference file
loci_file=paste0("~/pipeline_files/loci_files/loci_", build, ".RData")
load(loci_file)
rm(loci_file)
num_samples_pass <-  0

if(neo == TRUE) {
    # load neos
    df_SNVs <- data.table(
      read_csv("~/rdeocd/txMeth/neo_list_ext_new.csv"), stringsAsFactors = FALSE)
    df_SNVs <- df_SNVs[, .(Patient, Region, chr=chrom, start,
                           end, ref, var, mcn, vaf)]
  } else {
    # load all tx100 exome SNV calls
    df_SNVs <- data.table(read.delim(
      file="~/rdeocd/txMeth/tx_SNV_calls.txt", stringsAsFactors = FALSE))
    df_SNVs <- df_SNVs[GL_VAF < 0.1, ]
    df_SNVs <- df_SNVs[, .(SampleID, chr, start, end = stop, ref, var, MutCPN, RegionSum)]
    df_SNVs <- df_SNVs[, chr := as.character(chr)]
    df_SNVs <- df_SNVs[(chr != "24"), ]
    df_SNVs[chr == "23", chr := "X"]
    df_SNVs <- df_SNVs[end - start == 0]
}

# load patient_info
patient_info <- as.data.table(read.delim(
  paste0(parent_dir,
    "tx100clinicalinfo_RRBS_published_status_immune_clusters",
    "_purity_ploidy_CpG_cov.txt"),
  stringsAsFactors = FALSE))
patient_info <- patient_info[!grepl("LN", Sample), ]

# remove lymph nodes
patient_ids <- unique(patient_info$Patient)

for(patient_id in patient_ids[patient_ids != "LTX164"]) {
  path_patient <- paste0("~/rdeocd/txMeth/", patient_id, "/")
  path <- path_patient 
  gender <- as.character(
    gender_list[match(patient_id, gender_list$SampleID), "Gender"]
    )
  
  # set wd
  path_output=file.path(paste0(path_patient, "SNVs"))
  dir.create(path_output, showWarnings = FALSE)
  setwd(path_output)
 
  df_SNVs2 <- NULL
  f_nm <- paste0(path_patient, "SNVs/", patient_id, ".WGS.SNV.txt")
  if(file.exists(f_nm)&neo == FALSE){
    # load WGS SNVs
    df_SNVs2 <- data.table(read_tsv(f_nm, col_names = TRUE,
                                   col_types = cols(.default = "c")))
    # remove germline SNVs
    cols <- colnames(df_SNVs2)[grepl(".VAF", colnames(df_SNVs2))&grepl("GL", colnames(df_SNVs2))]
    idx <- which(as.numeric(data.frame(df_SNVs2, stringsAsFactors = FALSE)[, cols])>0)
    if(length(idx)>0){df_SNVs2 <- df_SNVs2[!c(1:nrow(df_SNVs2))%in%idx,]}
    rm(idx, cols)
    
    # keep relevant columns
    cols <- c("chr", "start", "ref", "var", 
              colnames(df_SNVs2)[grepl("[0-9].VAF$", colnames(df_SNVs2))&!grepl("GL", colnames(df_SNVs2))], 
              "is_SNV", "is_blacklist","Use.For.Plots")
    df_SNVs2 <- df_SNVs2[, .SD, .SDcols=cols]
    
    # filter SNVs
    df_SNVs2 <- df_SNVs2[Use.For.Plots==TRUE & is_blacklist==FALSE & is_SNV==TRUE, ]
    df_SNVs2[, c("is_SNV", "is_blacklist","Use.For.Plots") := NULL]
    
    # format
    colnames(df_SNVs2)[c(1,2)] <- c("chrom", "start")
    df_SNVs2$end <- df_SNVs2$start
  }
  rm(f_nm)
  
  # load heterozygous SNPs
  load(paste0("~/rdeocd/txMeth/",patient_id,"/Copy_number/N/",patient_id,".N.SNPs.RData"))
  df_het_SNPs <- data.table(df)[type == "Heterozygous", ]
  rm(df)

  # make heterozygous SNPs Granges
  colnames(df_het_SNPs)[c(1,3)] <- c("chrom", "start")
  df_het_SNPs$end <- df_het_SNPs$start
  gr_SNPs <- makeGRangesFromDataFrame(df_het_SNPs, keep.extra.columns = TRUE)
  rm(df_het_SNPs)
  
  if(neo == TRUE){
  #get patient SNVs
  df_SNVs_patient <- df_SNVs[Patient==patient_id, ]
  } else {
  #get patient SNVs
  CRUKID <- unique(patient_info[Patient == patient_id, CRUKID])
  df_SNVs_patient <- df_SNVs[SampleID==CRUKID, ]
  df_SNVs_patient <- df_SNVs_patient[, SampleID := NULL]
  rm(CRUKID)
  }
  
  #get sample_ids
  sample_ids <- unique(patient_info[Patient==patient_id, Sample])
  sample_ids <- sample_ids[!grepl("LN", sample_ids)]
  
  for(sample_id in sample_ids){
     dir.create(paste0("~/rdeocd/txMeth/",patient_id,"/SNVs/",sample_id))
      # # clean-up previous versions
      # f_nm <- paste0("~/rdeocd/txMeth/",patient_id,"/SNVs/", 
      #                sample_id,"/allele_specific_methylation_rates_",sample_id,".RData")
      # if(file.exists(f_nm)){
      #   file.remove(f_nm)
      # } ; rm(f_nm)
                 
      # sample pertaining to this loop iteration 
      print(paste(patient_id,sample_id,sep = "-"))
      
    if(neo == FALSE){
      #get sample-wise exome SNVs
      df_SNVs_sample <- df_SNVs_patient[grepl(sample_id, MutCPN),]
      df_SNVs_sample <- df_SNVs_sample[!is.na(MutCPN)]
    } else {
      df_SNVs_sample <- df_SNVs_patient[sample_id==Region,]
    }
    
  if(nrow(df_SNVs_sample)>1){
    if(neo == FALSE){
      # get SNV copy number per sample
      tmp <- data.table(str_split(df_SNVs_sample$MutCPN, ";", simplify = T))
      tmp$idx <- tmp[,apply(.SD, 1, function (x) grep(sample_id, x)), by=.(1:nrow(tmp))][,2]
      df_SNVs_sample[,mut_copy_number_string := tmp[,.SD[[idx]], by=.(1:nrow(tmp))][,2]]
      df_SNVs_sample[,mcn := as.numeric(as.character(str_split(mut_copy_number_string, ":", simplify = T)[,2]))]
      df_SNVs_sample[,mut_copy_number_string := NULL]
      rm(tmp)
      
      # get SNV vaf per sample
      tmp <- df_SNVs_sample[, tstrsplit(RegionSum,";", fixed=TRUE)]
      tmp$idx <- grep(sample_id,tmp[1])
      mut_vaf_string <- tmp[,.SD[[idx]], by=.(1:nrow(tmp))]
      colnames(mut_vaf_string)[2] <- "vaf_string"
      mut_vaf_string2 <- mut_vaf_string[, tstrsplit(vaf_string,":", fixed=TRUE)]
      vaf_num <- sapply(mut_vaf_string2$V2, function(x) eval(parse(text=x)))
      df_SNVs_sample[,vaf := as.numeric(vaf_num)]
      rm(tmp, mut_vaf_string, mut_vaf_string2, vaf_num)
      
      # rm SNVs with coverage = 0 
      df_SNVs_sample <- df_SNVs_sample[mcn>0 ,]
      df_SNVs_sample[,MutCPN := NULL]
    }
      # correct chrom anno
      df_SNVs_sample[, chrom := ifelse(substr(chr, 1, 1)!="c", paste0("chr", chr), chr)]
      df_SNVs_sample[, chr := NULL]
      
      # make GRanges
      gr_SNVs <- makeGRangesFromDataFrame(df_SNVs_sample, keep.extra.columns = TRUE)
      seqlevelsStyle(gr_SNVs) <- "UCSC"
      rm(df_SNVs_sample)
      
      if(!is.null(df_SNVs2)){
        #find sample_ids with WGS
        cols <- colnames(df_SNVs2)[grepl("[0-9].VAF", colnames(df_SNVs2))]
        sample_ids2 <- unlist(str_split(cols, "_"))[seq(3, length(cols)*3, 3)]
        sample_ids2 <- unlist(str_split(sample_ids2, "\\."))
        sample_ids2 <- sample_ids2[sample_ids2!="VAF"]
        
        if(sample_id%in%sample_ids2){
        # get relevant WGS VAF
        cols_sample <- colnames(df_SNVs2)[grepl(sample_id, colnames(df_SNVs2))&grepl("VAF", colnames(df_SNVs2))]
        cols <- c("chrom", "start", "end", "ref", "var", cols_sample)
        tmp <- df_SNVs2[, .SD, .SDcols=cols]
        colnames(tmp)[which(colnames(tmp)==cols_sample)] <- "vaf"
        
        # make Granges
        gr_SNVs2 <- makeGRangesFromDataFrame(tmp, keep.extra.columns = TRUE)
        seqlevelsStyle(gr_SNVs2) <- "UCSC"
        gr_SNVs2$mcn <- 0
        rm(tmp, cols_sample, cols)
        
        # make vaf numeric and fraction
        gr_SNVs2$vaf <- as.numeric(gr_SNVs2$vaf)/100
        
        # merge exome and WGS SNVs
        gr_SNVs <- sort(unlist(GRangesList(gr_SNVs, gr_SNVs2)))
        rm(gr_SNVs2)
        } ; rm(sample_ids2)
      } 
  
      # load copy number
      f_nm <- paste0(path_patient, "Copy_number/", sample_id,"/re_segmented/", patient_id,".", sample_id, ".ascat.output.resegmented.RData")
      if(file.exists(f_nm)){
        load(f_nm)
      } else {
        load(paste0(path_patient, "Copy_number/", sample_id,"/raw/", patient_id,".", sample_id, ".ascat.output.RData"))
      }
      
      #load copy number
      seg <- ascat.output$segments[,2:ncol(ascat.output$segments)]
      
      #find n+0 segments
      if(LOH_analysis==TRUE){
        seg <- seg[as.numeric(seg$nMinor)==0,]
      }
      #extract purity
      rho <- as.numeric(ascat.output$aberrantcellfraction)
      rm(ascat.output)
      
      #make Granges
      colnames(seg)[1:3] <- c("chrom", "start", "end")
      seg$chrom <- ifelse(substr(seg$chrom, 1, 1)!="c", paste0("chr", seg$chrom), seg$chrom)
      if(gender=="XY"){seg <- seg[!as.character(seg$chrom)%in%c("chrX", "chrY"),]}
      gr_seg <- makeGRangesFromDataFrame(seg, keep.extra.columns = TRUE)
      seqlevelsStyle(gr_seg) <- "UCSC"
      
      # annotate SNV copy number based on RRBS
      ov <- suppressWarnings(mergeByOverlaps(gr_SNVs, gr_seg))
      #x <- as.numeric(ov[, cols_sample]) >= (0.75 * ((ov$nMajor*rho)/(ov$nMajor*rho+2*(1-rho))))
      rm(gr_seg)
      
      # combine VAF and copy number info to get MCN where missing (i.e. WGS samples)
      df_SNV_loci <- data.table(data.frame(ov$gr_SNVs, nMinor=ov$nMinor, nMajor=ov$nMajor, stringsAsFactors = FALSE))
      df_SNV_loci[, nTotal := nMinor+nMajor]
      df_SNV_loci[df_SNV_loci$mcn==0|is.na(df_SNV_loci$mcn), mcn:= (1/rho)*vaf*nTotal*rho+2*(1-rho)]
      rm(ov)
      #plot(density(df_SNV_loci$mcn), xlim=c(0,4))
      
      # filter so mutation is on all copies of retained allele in tumour LOH regions 
      if(LOH_analysis==TRUE){
      df_SNV_loci <- df_SNV_loci[nMajor-0.3 <= mcn & nMajor+0.3 > mcn,]
      }
      # filter so mutation is likely clonal
      df_SNV_loci <- df_SNV_loci[vaf>0.1,]
      
    if(nrow(df_SNV_loci)>1){
        #x <- df_SNV_loci$VAF >= (0.75 * ((df_SNV_loci$nMajor*rho)/(df_SNV_loci$nMajor*rho+2*(1-rho))))
        #plot(density((0.75 * ((df_SNV_loci$nMajor*rho)/(df_SNV_loci$nMajor*rho+2*(1-rho))))), xlim=c(0,1))
     
      # make Granges for scanBAM
      gr_SNV_loci <-  makeGRangesFromDataFrame(df_SNV_loci, keep.extra.columns = TRUE)
      seqlevelsStyle(gr_SNV_loci) <- "UCSC"
      gr_SNV_loci$alt <- gr_SNV_loci$var
      gr_SNV_loci$POS <- start(gr_SNV_loci) 
      mcols(gr_SNV_loci) <- mcols(gr_SNV_loci)[, c("POS", "ref", "alt")] 
      gr_SNV_loci <- sort(gr_SNV_loci) 
 
      # get bam_file
      bam_file=paste0(bam_dir,patient_id,"/RRBS/Nugen/",patient_id,"_",sample_id,"_n6dupsRemoved_sorted.bam")
      if(patient_id == "LTX064") {bam_file=paste0(bam_dir,patient_id,"/RRBS/Nugen/",patient_id,"_",sample_id,"_n6dupsRemoved_firstmateonly_sorted_reheader.bam")}
      if((patient_id == "LTX001" & sample_id =="N")|(patient_id == "LTX180" & sample_id =="R2")|(patient_id == "LTX225" & sample_id =="R2")) {
        bam_file=paste0(bam_dir,patient_id,"/RRBS/Nugen/",patient_id,"_",sample_id,"_top_n6dupsRemoved_sorted.bam")
      }    
      
      # Load data
      if(sample_id=="N"){
        load(paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/gr_normal_m.RData"))
        gr <- gr_normal_m
        seqlevelsStyle(gr) <- "UCSC"
        rm(gr_normal_m)
      }
      if(sample_id!="N"){
        load(paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/gr_tumour_m.RData"))
        gr <- gr_tumour_m
        seqlevelsStyle(gr) <- "UCSC"
        rm(gr_tumour_m)
      }
      
      # Get CpGs which could be within 1 read length to neareast SNV
      #distances <- data.frame(mcols(distanceToNearest(gr, gr_SNV_loci)))[,1]
      #gr <- gr[distances<=110]
      gr_tmp <- suppressWarnings(extend(x = gr_SNV_loci, 110, 110))
      gr <- gr[gr%over%gr_tmp] 
      rm(gr_tmp)
      
  if(length(gr)>0){
      # Keep SNVs which could be within 1 read length to neareast CpG
      gr_tmp <- suppressWarnings(extend(x = gr, 110, 110))
      gr_SNV_loci <- gr_SNV_loci[gr_SNV_loci%over%gr_tmp] 
      rm(gr_tmp)
      
      # remove CpGs overlapping with SNVs
      gr <- gr[!gr%over%gr_SNVs] 
      gr <- gr[!duplicated(gr)]
      
      # remove CpGs overlapping with het SNPs
      gr <- gr[!gr%over%gr_SNPs]
      
      # Set build variables
      if(build%in%c("hg19","hg38")){UCSC=TRUE}
      if(build%in%c("GRCH37","GRCH38")){UCSC=FALSE}
      if(build=="GRCH37"){build="hg19"} # set build to to assembly version disregarging USCS/Ensembl
      if(build=="GRCH38"){build="hg38"}
      cat(paste("ScanBam pileup with build ", build, sep = " "), "\n", sep = "")
      
      # Only CpGs near SNVs but not at SNVs
      loci_subset <- all_loci[all_loci%in%gr]
      loci_subset <- loci_subset[!loci_subset%over%gr_SNVs]
      rm(gr)
      
      # Only keep CpGs outside het SNPs
      loci_subset <- loci_subset[!loci_subset%over%gr_SNPs]
      loci_subset <- loci_subset[width(loci_subset)>1]
      
      # remove SNPs info at hom SNPs
      loci_subset$POS <- NA 
      loci_subset$ref <- NA 
      loci_subset$alt <- NA 
      
  if(length(loci_subset)>=3&length(gr_SNV_loci)>0){
      # add SNVs to loci_file
      loci_subset <- sort(unlist(GRangesList(loci_subset, gr_SNV_loci)))
      names(loci_subset) <- 1:length(loci_subset)
      
  if(length(loci_subset)>=3){
      # Set what argument to speed up scanBam
      whats <- c("qname", "rname", "strand", "pos", "qwidth", "mapq", "seq", "qual") # scanBamWhat()
      idxFile <- paste(bam_file, ".bai", sep = "")
      
      # Set info that we want to obtain for each position
      bamFile <- BamFile(bam_file)
      params <- Rsamtools::ScanBamParam(which=loci_subset, what=whats)
      
      # Run ScanBam 
      #bam <- Rsamtools::scanBam(file=bam_file, index=idxFile, param=params)
      bam <- Rsamtools::scanBam(file=bamFile, param=params)
      rm(params)
      
      # Remove entries with zero read
      bam_filtered <- bam[which(sapply(bam, function(x) length(x[["qname"]])) != 0)]
      rm(bam)
      
      # Unlist bam
      df <- data.table(dplyr::bind_rows(lapply(bam_filtered, as.data.frame.list, stringsAsFactors=FALSE)))
      # duplicates created when unlisting and due to overlapping reads, to be removed in later lines
      rm(bam_filtered)
      
      # Format pileup df for transformation into GRanges Object
      df$end <- df$pos + (df$qwidth - 1)
      colnames(df)[c(which(colnames(df) == "rname"), 
                     which(colnames(df) == "pos"))] <- c("chrom", "start")
      rownames(df) <- c(1:nrow(df))
      
      cols <- c("qname", "chrom", "strand", "start", "end", "mapq", "seq", "qual")
      Sys.sleep(1)
      
      df_bam <- df[,cols,with=FALSE] 
      rm(df)
      
      # The pileup unlisting produces duplicate rows where more than one loci overlaps with the same molecule
      x <- do.call(paste, c(df_bam[, cols, with=FALSE], sep="-")) 
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
     
      cols<-c("qname","chrom","POS","width","start","ref","alt")
      x <- do.call(paste, c(df_pileup[,cols,with=FALSE], sep="-")) 
      df_pileup_dedup <- df_pileup[!duplicated(x),]
      
      # clean-up
      rm(df_pileup, cols, x, overlaps)
      
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
      tmp[, ID := paste(chrom, POS, ref, alt, sep="_")]
      tmp2 <- paste(seqnames(gr_SNV_loci), gr_SNV_loci$POS, gr_SNV_loci$ref, gr_SNV_loci$alt, sep="_")
      # tmp[, table(ID%in%as.character(tmp2))]
      tmp <- tmp[ID%in%as.character(tmp2),]
      rm(tmp2)
      
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
      tmp <- tmp[,.("ref_counts" = sum(ref_counts, na.rm = TRUE), 
                    "alt_counts" = sum(alt_counts, na.rm = TRUE)),
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
                               POS_genotype = paste0(POS,collapse = ","), 
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
      
     if(sum(!is.na(df_pileup_alleles$genotype)&is.na(df_pileup_alleles$POS)) > 0){
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
        
    if(nrow(flag)>0){
        # Keep position with clear ref / alt alleles
        df_pileup_alleles <- merge(df_pileup_alleles, flag[,.(chrom, start, end, width, POS, ref, alt)], 
                                   by = c("chrom", "POS", "start", "end", "width","ref", "alt"))
        
        # Calculate total ref or alt count for each position
        df_summary <- df_pileup_alleles[,.("Af" = sum(Af, na.rm = TRUE), "Ar" = sum(Ar, na.rm = TRUE),"Cf" = sum(Cf, na.rm = TRUE), "Cr" = sum(Cr, na.rm = TRUE),
                                           "Gf" = sum(Gf, na.rm = TRUE), "Gr" = sum(Gr, na.rm = TRUE),"Tf" = sum(Tf, na.rm = TRUE), "Tr" = sum(Tr, na.rm = TRUE),
                                           "CGf" = sum(CGf, na.rm = TRUE), "CGr" = sum(CGr, na.rm = TRUE),"TGf" = sum(TGf, na.rm = TRUE), "CAr" = sum(CAr, na.rm = TRUE),
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
        
        # Where there is insufficient evidence that a CG motif is present or what its methylation status is, set M and UM to NA
        df_summary[total_counts_m <= 3 & width >= 2, c("M","UM") := NA ]
        
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
        flag <- df_summary[, (ifelse(other_counts <= 0.05*total_depth, FALSE, TRUE))]
        if(sum(flag)>0) {
          warning(paste0("Unexpected nucleotides represent > 5% of the alleles at ", sum(flag), 
                         " positions ( ", round(sum(flag)*100/length(flag), digits = 2),
                         " %) - positions removed.", sep = ""))
          df_summary<- df_summary[!flag,]
        };rm(flag)
        
        df <- df_summary[,c("chrom","start","end","width","POS","ref","alt","alt_counts","ref_counts",
                            "total_counts","total_depth","other_counts","all_counts","M","UM","total_counts_m",
                            "Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr","CCGG",
                            "genotype", "POS_genotype", "SNP"), with = F]
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
        if (getOption("scipen") == 0) {options(scipen = 999)}
        
        # Ensure there are no NAs before saving
        flag <- df[,!is.na(CHR)]
        if(sum(flag==FALSE)>0){df <- df[flag,]};rm(flag)
        
        # Select columns for saving and put in a logical order
        df <- df[,c("CHR","chrom","start","end","width","POS","ref","alt",
                    "alt_counts","ref_counts","total_counts","BAF","total_depth",
                    "other_counts","all_counts","M","UM","total_counts_m","m",
                    "Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr","CCGG",
                    "genotype", "POS_genotype", "SNP"), 
                 with=F]
        
        cat(paste("ASM pileup completed", build, sep = " "), "\n", sep = "")
        
        # keep rows with both alleles remaining after QC
        flag <- df[,(num_genotypes = length(unique(na.omit(genotype)))), 
                   keyby=c("CHR","chrom","start","end","width","POS","ref","alt")]
        flag <- flag[V1==2,]
        df <- merge(df,flag,by=c("CHR","chrom","start","end","width","POS","ref","alt"))
        
        # check matched normal exists
        fls = paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/N/gr_normal_m.RData")
        if(!file.exists(fls)){
            warnings("no data for matched normal")
            stopifnot(file.exists(fls))
        } else {
           load(fls)
        }
          
        if(file.exists(fls)){
          # annotate CN
          seg <- data.table(seg)
          colnames(seg)[1:3] = c("CHR", "start", "end")
          setkeyv(seg, c("CHR", "start", "end"))
          setkeyv(df,  c("CHR", "start", "end"))
          tmp <- foverlaps(df, seg, type = "within") 
          df$nA =  tmp$nMajor[match(paste(df$CHR, df$POS_genotype, sep="_"), 
                                    paste(tmp$CHR, tmp$POS_genotype, sep="_"))]
          df$nB = tmp$nMinor[match(paste(df$CHR, df$POS_genotype, sep="_"), 
                                   paste(tmp$CHR, tmp$POS_genotype, sep="_"))]
          df_tumour <- df
          df_tumour <- df_tumour[!is.na(total_counts_m),]
          rm(df, seg)
          
          if(nrow(df_tumour)>=2){
            num_samples_pass <-  num_samples_pass +1
            # load matched normal
            seqlevelsStyle(gr_normal_m) <- "Ensembl"
            df <- data.table(data.frame(gr_normal_m))
            colnames(df)[1] <- "chrom"
            df_normal <- df[,.(chrom, start, end, M_n, UM_n, total_counts_m_n = M_n+UM_n)]
            rm(df, gr_normal_m)
                
            # join tumour and normal info
            df_tumour <- merge(df_tumour, df_normal, by = c("chrom", "start", "end"), all.x=T) 
            rm(df_normal)
            rho <- as.numeric(rho) 
            
            # get allele-specific methylation
            df_phased_methylation <- df_tumour[order(chrom,start,end,width,POS,genotype), 
                                               .(CN=nA+nB,genotype_1="ref",genotype_2="alt", num_genotypes=length(genotype),
                                                POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                SNV=SNP[genotype=="ref"&SNP%in%SNP[genotype=="alt"]],
                                                counts_1=total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                counts_2=total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                                M_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                                UM_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                total_counts_m_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                M_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                                UM_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                                total_counts_m_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                                M_n = M_n[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                UM_n = UM_n[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                                total_counts_m_n = total_counts_m_n[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]]), 
                                               keyby=c("CHR","chrom","start","end","width","POS","ref","alt", "nA", "nB")]  
            #rm(df_tumour)
            
            # get phased tumour and bulk normal methylation rates    
            df_phased_methylation[, m_1 := as.numeric(M_1)/as.numeric(total_counts_m_1)] 
            df_phased_methylation[, m_2 := as.numeric(M_2)/as.numeric(total_counts_m_2)]
            df_phased_methylation[, m_n := as.numeric(M_n)/as.numeric(total_counts_m_n)] 
            
            # observed read counts ratio with and without SNV
            df_phased_methylation[, SNV_VAF_methylation := as.numeric(total_counts_m_2)/
                                    (as.numeric(total_counts_m_2)+as.numeric(total_counts_m_1))] 
                
            # methylation mcn
            df_phased_methylation[, SNV_mcn_methylation := (1/rho)*SNV_VAF_methylation*((nA+nB)*rho+2*(1-rho))] 

            # # expected read ratio ratio with and without SNV from WGS
            # cols_sample <- colnames(df_SNVs)[grepl(sample_id, colnames(df_SNVs))&grepl("VAF", colnames(df_SNVs))]
            # df_SNVs_tmp <- df_SNVs[, .SD, .SDcols=c("chrom", "start", cols_sample)]
            # df_SNVs_tmp[, SNV_VAF := as.numeric(unlist(df_SNVs_tmp[, cols_sample, with=FALSE]))/100] 
            
            # expected read ratio ratio with and without SNV from WGS
            df_SNVs_tmp <- df_SNV_loci[, .SD, .SDcols=c("seqnames", "start", "vaf")]
            df_SNVs_tmp[, SNV_VAF := vaf]
            
            # assign to df
            df_phased_methylation$SNV_VAF_WGS <- df_SNVs_tmp[match( 
              paste(df_phased_methylation$CHR, df_phased_methylation$POS_genotype, sep="_"), 
              paste(seqnames, start, sep="_")),]$SNV_VAF
            
            # expected read ratio ratio with and without SNV from WGS
            df_SNVs_tmp <- df_SNV_loci[, .SD, .SDcols=c("seqnames", "start", "mcn")]
            df_SNVs_tmp[, SNV_mcn := mcn]
            
            # assign to df
            df_phased_methylation$SNV_mcn <- df_SNVs_tmp[match( 
              paste(df_phased_methylation$CHR, df_phased_methylation$POS_genotype, sep="_"), 
              paste(seqnames, start, sep="_")),]$SNV_mcn
            
            # # expected read ratio ratio with and without SNV from model
            # df_phased_methylation[, SNV_VAF_model := ((nA+nB)*rho)/((nA+nB)*rho+2*(1-rho))]
            # rm(df_SNVs_tmp)
            
            # carry out tumour-normal deconvolution
            df_phased_methylation <- compute_tumour_methylome(df_phased_methylation, p=rho, min_cov_t = 3, gender = gender)
            cat(paste("ASM deconvolution completed", build, sep = " "), "\n", sep = "")
            
            # at LOH regions tumour-normal decon is based on allele phasing only
            df_phased_methylation[nB==0, m_t_1 := NA]
            df_phased_methylation[nB==0, m_t_2 := m_2]
            
            # get potentially imprinted loci
            df_phased_methylation[,m_diff_t12 := abs(m_t_1-m_t_2)]
            df_phased_methylation[,m_diff_b12 := abs(m_1-m_2)]
            
            # get methylation difference
            df_phased_methylation[,m_diff_tn1 := abs(m_n-m_t_1)]
            df_phased_methylation[,m_diff_tn2 := abs(m_n-m_t_2)]
            df_phased_methylation[,m_diff_bn1 := abs(m_n-m_1)]
            df_phased_methylation[,m_diff_bn2 := abs(m_n-m_2)]
            df_phased_methylation[,m_diff_tn := abs(m_n-m_t)]
            df_phased_methylation[,m_diff_bn := abs(m_n-m_b)]
             
            # remove NAs
            df_phased_methylation <- df_phased_methylation[!is.na(m_t)&!is.na(m_t_2),]
            df_phased_methylation <- df_phased_methylation[!duplicated(df_phased_methylation),]
            
            #  save data table
            dir_path <- file.path(path_patient, "SNVs", sample_id)
            suppressWarnings(dir.create(dir_path, showWarnings = FALSE))
            if(neo==FALSE){
            save(df_phased_methylation, file=paste0(dir_path,"/allele_specific_methylation_rates_",sample_id,".RData"))
              f_nm <- paste0(dir_path,"/ASM_SNVs.pdf")
            } else {
            save(df_phased_methylation, file=paste0(dir_path,"/neo_allele_specific_methylation_rates_",sample_id,".RData"))
              f_nm <- paste0(dir_path,"/ASM_neos.pdf")
            }
            max_plus <- max(na.omit(c(df_phased_methylation$m_t, df_phased_methylation$m_t_2)))
            max_plus <- ifelse(max_plus < 1.45, 1.45, max_plus )
            max_minus <- max(na.omit(c(-df_phased_methylation$m_t, -df_phased_methylation$m_t_2)))
            max_minus <- ifelse(max_minus < 0.5, 0.5, max_minus)
            
            if(file.exists(f_nm)){file.remove(f_nm)}
            if(nrow(df_phased_methylation)>10){
            p1 <- ggplot(df_phased_methylation, aes(x=m_t_2, y=m_t)) + 
              geom_point(size=2, shape=16, alpha=0.5, col="navyblue") + theme_classic() + 
              geom_smooth(method="lm") +
              ylim(-max_minus,max_plus) + xlim(-max_minus,max_plus) + 
              geom_text(data=df_phased_methylation[1,], aes(x = -0.5, y = 1.2), size=2, 
                        label=paste0("n = ", nrow(df_phased_methylation)), vjust=0, hjust=0)
            p1 <- p1 + 
              geom_text(data=df_phased_methylation[1,], aes(x = -0.5, y = 1.4), size=2, 
                        label = lm_eqn(x=df_phased_methylation$m_t_2, y=df_phased_methylation$m_t, 
                                       df=df_phased_methylation), 
                                 parse = TRUE, vjust=0, hjust=0)
            tmp <- data.frame(flag=rep(c("Tmut_v_TWT", "Tmut_v_NWT", "TWT_v_NWT"), each=nrow(df_phased_methylation)),
                              m_diff = c(df_phased_methylation$m_diff_t12, df_phased_methylation$m_diff_tn2,
                                         df_phased_methylation$m_diff_tn1))
            tmp <- tmp[!is.na(tmp$m_diff),]
            p2 <- ggplot(tmp, aes(x=m_diff, y=..density.., col=flag)) + 
              geom_density() + theme_classic() + 
              coord_cartesian(xlim=c(-0.4,1.4)) 
            
            # set grid
            gglist1 <- list(p1) 
            gglist2 <- list(p2) 
            my_layout <- cbind(1,1,1,1,2,2,2,2,2)
            gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
            ggsave(filename =  f_nm, plot=gd, device = "pdf", units = "in", width = 7, height =2.5)
            rm(p1, p2, g2, gglist1, gglist2)
          } ; rm(max_plus, max_minus, dir_path, f_nm)
          }
        } 
        }
       }
      }
     }
    }
   }
  }
 } 
}
rm(df_tumour, all_loci, df_SNVs, df_SNVs2, df_SNV_loci, gr_SNV_loci, gr_SNVs, gr_SNPs)

list_patients <- list()
for(patient_id in patient_ids[patient_ids!="LTX164"]) {
  path_patient <- paste0("~/rdeocd/txMeth/",patient_id,"/")
  path = path_patient 
  gender=as.character(gender_list[match(patient_id,gender_list$SampleID),"Gender"])
  
  # set wd
  path_output=file.path(path_patient, "SNVs")
  setwd(path_output)
  
  # # get sample ids
  # sample_ids <- list.files(file.path(path_patient, "WGS")) 
  # sample_ids <- sample_ids[grep("^R",sample_ids)]
  
  #get sample_ids
  sample_ids <- patient_info[Patient==patient_id, unique(Sample)]
  
  list_samples <- list()
  for (sample_id in sample_ids){
    f_nm = file.path(path_patient, "SNVs", sample_id, paste0("allele_specific_methylation_rates_", sample_id,".RData"))
    #f_nm = file.path(path_patient, "SNVs", sample_id, paste0("neo_allele_specific_methylation_rates_", sample_id,".RData"))
    if(file.exists(f_nm)){
      load(f_nm)
      df_phased_methylation[, LTX_ID := paste(patient_id, sample_id, sep="_")]
      list_samples[[length(list_samples)+1]] <- df_phased_methylation
    } 
  }
  if(length(list_samples)>0){
  list_patients[[length(list_patients)+1]] <- rbindlist(list_samples)
  }
  rm(list_samples)
}

# bind list of data.frame with each sample's informations
all_df_phased_methylation <- rbindlist(list_patients, fill=TRUE)
rm(list_patients)
save(all_df_phased_methylation, file="~/rdeocd/txMeth/SNV_DMPs.RData")

# remove NAs and dups and SNVs not on all copies
all_phased_cpgs <- nrow(all_df_phased_methylation)
all_snvs <- nrow(
  unique(all_df_phased_methylation[, .SD, .SDcols = c("CHR", "POS_genotype", "LTX_ID")])
  )
df_phased_methylation <- # subset to regions of LOH
  all_df_phased_methylation[!is.na(m_t) & !is.na(m_t_2) & nB == 0,] 
df_phased_methylation <- # subset to regions where the SNVs is on all tumour copies
  df_phased_methylation[SNV_mcn_methylation >= 0.5 &
  abs(SNV_mcn_methylation - nA) <= 0.5]
nrow(
  unique(df_phased_methylation[, .SD, .SDcols = c("CHR", "POS_genotype", "LTX_ID")])
  )
df_phased_methylation <- df_phased_methylation[!duplicated(df_phased_methylation)]

# ge the number of unique mutations used

# plot corrected values
df_phased_methylation[, m_t_corr := ifelse(m_t < 0, 0, ifelse(m_t > 1, 1, m_t))]
df_phased_methylation[, m_t_2_corr := ifelse(m_t_2 < 0, 0, ifelse(m_t_2 > 1, 1, m_t_2))]

# add purity values
df_phased_methylation$rho <- 
  patient_info[match(df_phased_methylation$LTX_ID, LTX_ID), "rho_reseg"]
    
df_phased_methylation$SNV_VAF_WGS[df_phased_methylation$SNV_VAF_WGS==0] <- 
  (df_phased_methylation$SNV_mcn[df_phased_methylation$SNV_VAF_WGS==0]*
   df_phased_methylation$rho[df_phased_methylation$SNV_VAF_WGS==0])/
  (df_phased_methylation$nA[df_phased_methylation$SNV_VAF_WGS==0]*
   df_phased_methylation$rho[df_phased_methylation$SNV_VAF_WGS==0]+
   2*(1-df_phased_methylation$rho[df_phased_methylation$SNV_VAF_WGS==0]))

# set plot boundaries
xmax_plus <- df_phased_methylation[, max(SNV_mcn, na.rm=TRUE)]
xmax_minus <- df_phased_methylation[, min(SNV_mcn, na.rm=TRUE)]
ymax_plus <- df_phased_methylation[, max(SNV_mcn_methylation, na.rm=TRUE)]
ymax_minus <- df_phased_methylation[, min(SNV_mcn_methylation, na.rm=TRUE)]

# compare observed vs expected mutation copy number
p1 <- ggplot(df_phased_methylation, aes(x=SNV_mcn, y=SNV_mcn_methylation)) + 
  theme_classic() + xlab("Observed (methylation) - Expected (WGS/WES)") + 
  ylab("Observed (methylation) - Expected (WGS/WES)") +
  geom_point(size=1, shape=16, alpha=0.5, col="green3") +
  ylim(ymax_minus,ymax_plus) + xlim(xmax_minus,xmax_plus) + 
  geom_text(data=df_phased_methylation[1,], aes(x = xmax_minus+0.1, y = ymax_plus-0.1), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), 
                        "\ncor = ", round2(cor(df_phased_methylation[, SNV_mcn], 
                                               df_phased_methylation[, SNV_mcn_methylation]), 
                         digits = 2)), hjust=0)

# p1 <- p1 + geom_text(data=df_phased_methylation[1,], aes(x = xmax_minus+0.1, y = ymax_plus-0.2), size=2, 
#                      label = lm_eqn(x=df_phased_methylation[, (SNV_VAF_WGS-SNV_VAF_model)], 
#                                     y=df_phased_methylation[, (SNV_VAF_methylation-SNV_VAF_model)], 
#                                     df=df_phased_methylation), parse = TRUE, hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_mcn_RRBS_vs_GS_points.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

df_phased_methylation[, vaf_diff := abs(SNV_VAF_methylation - SNV_VAF_WGS)]
# compare observed vs expected vaf
p1 <- ggplot(df_phased_methylation, aes(x=SNV_VAF_WGS, y=SNV_VAF_methylation)) + 
  theme_classic() + xlab("WGS/WES VAF") + ylab("RRBS VAF") +
  #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
  geom_point(size=1, shape=16, alpha=0.25, col="black", fill="black") +
  #geom_smooth(method="lm", col="grey40") +
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="grey40") +
  ylim(0,1) + xlim(0,1) + scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") +
  # geom_hline(aes(yintercept=-0.25), col="grey60", size=0.2, linetype="dashed") +
  # geom_hline(aes(yintercept=0.25), col="grey60", size=0.2, linetype="dashed") + 
  geom_text(data=df_phased_methylation[1,], aes(x = 0.1, y = 0.8), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), "\ncor = ", round2(cor(
                  df_phased_methylation[, SNV_VAF_WGS], df_phased_methylation[, SNV_VAF_methylation]), 
            digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_vaf_RRBS_vs_GS_points.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

df <- data.frame(y = c(df_phased_methylation[,(SNV_mcn_methylation-SNV_mcn)]),
                 y2 = c(df_phased_methylation[,(SNV_VAF_methylation-SNV_VAF_WGS)]),
                 x = factor(c(rep(0, nrow(df_phased_methylation))), levels = c(0), ordered = TRUE), 
                 flag = c(rep("RRBS-GS", nrow(df_phased_methylation))))
# my_comparisons <- list(factor(c(0, 1), levels = c(0,1), ordered = TRUE))

p1 <- ggplot(df, aes(x=x, y=y, col=flag, fill=flag)) +  theme_classic() +  ylab("Observed - Expected")  +
  geom_boxplot(width=.75, alpha=0.1) +
  scale_colour_manual(name="", values=c("RRBS-GS"="darkorange")) +
  scale_fill_manual(name="", values=c("RRBS-GS"="darkorange")) + 
  #stat_compare_means(comparisons = my_comparisons, size = 2) + 
  scale_x_discrete(name="Platform", breaks = 0, labels = c("RRBS-GS"))
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_mcn_RRBS_vs_GS_box.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

p1 <- ggplot(df, aes(x=x, y=y2, col=flag, fill=flag)) +  theme_classic() +  
  ylab("vaf (observed - expected)")  +
  geom_boxplot(width=.75, alpha=0.1) +
  scale_colour_manual(name="", values=c("RRBS-GS"="darkorange")) +
  scale_fill_manual(name="", values=c("RRBS-GS"="darkorange")) + 
  #stat_compare_means(comparisons = my_comparisons, size = 2) 
  scale_x_discrete(name="Platform", breaks = 0, labels = c("RRBS-GS"))
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_vaf_RRBS_vs_GS_box.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

ggplot(data=df, aes(x=y, y=..density.., col=flag, fill=flag)) + theme_classic() + 
  xlab("Observed - Expected") +
  geom_density(alpha=0.1) + 
  scale_colour_manual(name="", values=c("RRBS-GS"="darkorange")) +
  scale_fill_manual(name="", values=c("RRBS-GS"="darkorange"))+ 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_mcn_RRBS_vs_GS_density.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

p1 <- ggplot(df, aes(x=y2, y=..density.., col=flag, fill=flag)) + theme_classic() + 
   xlab("vaf (observed - expected)") +
   geom_density(alpha=0.1) + 
   scale_colour_manual(name="", values=c("RRBS-GS"="darkorange")) +
   scale_fill_manual(name="", values=c("RRBS-GS"="darkorange"))
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_vaf_RRBS_vs_GS_density.pdf"), 
        plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

df_phased_methylation[,dist_SNV := abs(as.numeric(df_phased_methylation$POS_genotype)-df_phased_methylation$start)]
df_phased_methylation[,flag := dist_SNV<4 &width==4]
tmp <- df_phased_methylation[, .(CCGG = sum(flag)>0), by = .(chrom, POS_genotype)]
values <- tmp[match(paste0(df_phased_methylation$chrom, df_phased_methylation$POS_genotype), 
                     paste0(chrom, POS_genotype)), CCGG]
df_phased_methylation[,flag2 := values]
df <- data.frame(y2 = c(df_phased_methylation[flag2 ==FALSE,(SNV_VAF_methylation-SNV_VAF_WGS)], 
                        df_phased_methylation[flag2 ==TRUE,(SNV_VAF_methylation-SNV_VAF_WGS)]),
                 x = factor(c(rep(0, sum(df_phased_methylation$flag2 ==FALSE)),
                              rep(1, sum(df_phased_methylation$flag2==TRUE))), 
                              levels = c(0,1), ordered = TRUE), 
                 flag = c(rep("non-CCGG", sum(df_phased_methylation$flag2==FALSE)),
                          rep("CCGG", sum(df_phased_methylation$flag2==TRUE))))
p1 <- ggplot(df, aes(x=x, y=y2, col=flag, fill=flag)) +  theme_classic() +  
  ylab("vaf (observed - expected)")  +
  geom_boxplot(width=.75, alpha=0.1) +
  scale_colour_manual(name="", values=c("CCGG"="mediumpurple2","non-CCGG" = "grey50")) +
  scale_fill_manual(name="", values=c("CCGG"="mediumpurple2","non-CCGG" = "grey50")) + 
  stat_compare_means(size = 2) +
  scale_x_discrete(name="Platform", breaks = c(0,1), labels = c("non-CCGG","CCGG"))
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_vaf_RRBS_vs_GS_with_CCGG_context_box.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

# compare observed vs expected vaf - WGS
p1 <- ggplot(df_phased_methylation, aes(x=SNV_VAF_WGS, y=SNV_VAF_methylation)) + 
  theme_classic() + xlab("WGS/WES VAF") + ylab("RRBS VAF") +
  #geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) +
  geom_point(data=df_phased_methylation[flag2==FALSE], 
             aes(col=flag2, fill=flag2), size=1, shape=16, alpha=0.2) +
  geom_point(data=df_phased_methylation[flag2==TRUE], 
             aes(col=flag2, fill=flag2), size=1.5, shape=16, alpha=1) +
  #geom_smooth(method="lm", col="grey40") +
  scale_colour_manual(name="", values=c("FALSE"="grey50","TRUE"= "mediumpurple2"),
                      labels=c("FALSE"="non-CCGG", "TRUE"="CCGG")) +
  scale_fill_manual(name="", values=c("FALSE"="grey50","TRUE"="mediumpurple2"),
                    labels=c("FALSE"="non-CCGG", "TRUE"="CCGG")) + 
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="grey40") +
  ylim(0,1) + xlim(0,1) +
  # geom_hline(aes(yintercept=-0.25), col="grey60", size=0.2, linetype="dashed") +
  # geom_hline(aes(yintercept=0.25), col="grey60", size=0.2, linetype="dashed") + 
  geom_text(data=df_phased_methylation[1,], aes(x = 0.1, y = 0.8), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), 
                         "\ncor = ", round2(cor(
              df_phased_methylation[, SNV_VAF_WGS], 
              df_phased_methylation[, SNV_VAF_methylation]), 
              digits = 2),
              "\ncor without CCGG = ", round2(cor(
                df_phased_methylation[flag2==FALSE, SNV_VAF_WGS], 
                df_phased_methylation[flag2==FALSE, SNV_VAF_methylation]), 
              digits = 2),
              "\ncor at CCGG = ", round2(cor(
                df_phased_methylation[flag2==TRUE, SNV_VAF_WGS], 
                df_phased_methylation[flag2==TRUE, SNV_VAF_methylation]), 
              digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_vaf_RRBS_vs_GS_points.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

# set plot boundaries
max_plus <- max(na.omit(c(df_phased_methylation$m_t, df_phased_methylation$m_t_2)))
max_minus <- max(na.omit(c(-df_phased_methylation$m_t, -df_phased_methylation$m_t_2)))
 
df <- data.frame(y = c(df_phased_methylation[,(m_t_2)], df_phased_methylation[,(m_t)]),
                 x = factor(c(rep(0, nrow(df_phased_methylation)),rep(1, nrow(df_phased_methylation))), levels = c(0,1), ordered = TRUE), 
                 flag = c(rep("SNV", nrow(df_phased_methylation)),rep("CAMDAC", nrow(df_phased_methylation))))
my_comparisons <- list(factor(c(0, 1), levels = c(0,1), ordered = TRUE))

p1 <- ggplot(df, aes(x=x, y=y, col=flag, fill=flag)) +  theme_classic() +  
  ylab("pure tumour methylation rate")  +
  geom_boxplot(width=.75, alpha=0.1) +
  scale_colour_manual(name="", values=c("SNV"="grey30", "CAMDAC"="darkorange")) +
  scale_fill_manual(name="", values=c("SNV"="grey30", "CAMDAC"="darkorange")) + 
  stat_compare_means(comparisons = my_comparisons, size = 2) + 
  scale_x_discrete(name="Platform", breaks = seq(0, 1, 1), labels = c("SNV","CAMDAC"))
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_CAMDAC_box.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.5, height =2.8)

p1 <- ggplot(df_phased_methylation, aes(x=m_t_2, y=m_t)) + 
  geom_point(size=0.5, shape=16, alpha=0.25, col="navyblue") + theme_classic() + 
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
  #geom_smooth(method="lm", formula = y ~ 0 + x) +
  ylim(-max_minus,max_plus) + xlim(0,1) + 
  geom_text(data=df_phased_methylation[1,], aes(x = 0, y = 1.35), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), 
                         "\ncor = ", round2(cor(df_phased_methylation[, m_t], 
                                       df_phased_methylation[, m_t_2]), digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_CAMDAC_decon.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

# lm_eqn <- function(df, x, y){
#   m <- lm(y ~ 0 + x, df);
#   eq <- substitute(italic(y) == a %.% italic(x)*","~~italic(r)^2~"="~r2,
#                    list(a = format(unname(coef(m)[1]), digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
# p1 <- p1 + geom_text(data=df_phased_methylation[1,], aes(x = 0.3, y = 0.75), size=2,
#                      label = lm_eqn(x=c(df_phased_methylation$m_t),
#                                     y=c(df_phased_methylation$m_t_2),
#                                     df=df_phased_methylation), parse = TRUE, vjust=0, hjust=0) 

p1 <- ggplot(df_phased_methylation[SNV_VAF_WGS!=0, ], 
             aes(x=m_t-m_t_2, y=SNV_VAF_methylation-SNV_VAF_WGS)) +
  geom_point(size=1, shape=16, alpha=0.5, col="navyblue") + theme_classic() 

ggplot(df_phased_methylation, aes(x=m_2-m_t, y=..density..)) + 
  geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() + 
  geom_text(data=df_phased_methylation[1,], aes(x = -max_minus+0.5, y = 3), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), 
                         "\ncor = ", round2(cor(df_phased_methylation[, m_b], 
                                                df_phased_methylation[, m_2]), digits = 2)), hjust=0)

p1 <- ggplot(df_phased_methylation, aes(x=m_t_2-m_t, y=..density..)) + 
  geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() + 
  geom_text(data=df_phased_methylation[1,], aes(x = -max_minus+0.5, y = 3), size=2, 
            label=paste0("n = ", nrow(df_phased_methylation), 
                         "\ncor = ", round2(cor(df_phased_methylation[, m_t], 
                                                df_phased_methylation[, m_t_2]), digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_CAMDAC_decon_density.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

tmp <- df_phased_methylation
p <- ggplot(tmp, aes(x=m_2, y=m_b)) + theme_classic() + 
  geom_point(size=1, shape=16, alpha=0.25, col="black", fill="black") +
  ylim(0,1) + xlim(0,1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="black") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
            label=paste0("n = ", nrow(tmp), "\ncor = ", 
                         round2(cor(tmp[, m_b], tmp[, m_2]),  digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_bulk.pdf"), 
       plot=p, device = "pdf", units = "in", width = 3.4, height =2.8)

p1 <- ggplot(tmp, aes(x=m_t_2_corr, y=m_t_corr)) + theme_classic() + 
  geom_point(size=1, shape=16, alpha=0.25, col="navyblue") +
  ylim(0,1) + xlim(0,1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="black") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
     label=paste0("n = ", nrow(tmp), "\ncor = ", 
                  round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]),  digits = 2)), hjust=0)
# p1 <- p1 + geom_text(data=tmp[1,], aes(x = 0, y = 0.75), size=2, 
#                      label = lm_eqn(x=tmp$m_t_2, y=tmp$m_t, 
#                                     df=tmp), parse = TRUE, vjust=0, hjust=0) 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_CAMDAC_decon_capped.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

tmp <- df_phased_methylation[abs(m_diff_tn2) > 0.2]
p <- ggplot(tmp, aes(x=m_2, y=m_b)) + theme_classic() + 
  geom_point(size=1, shape=16, alpha=0.25, col="black") +
  ylim(0,1) + xlim(0,1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="black") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
            label=paste0("n = ", nrow(tmp), "\ncor = ", 
                         round2(cor(tmp[, m_b], tmp[, m_2]),  digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_bulk_tn2_DMP.pdf"), 
       plot = p, device = "pdf", units = "in", width = 3.4, height = 2.8)

p1 <- ggplot(tmp, aes(x=m_t_2_corr, y=m_t_corr)) + theme_classic() + 
  geom_point(size=1, shape=16, alpha=0.25, col="navyblue") + # aes(col=vaf_diff)
  ylim(0,1) + xlim(0,1) + 
  #scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") + 
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
            label=paste0("n = ", nrow(tmp), "\ncor = ", 
                         round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]),  digits = 2)), hjust=0)
# p1 <- p1 + geom_text(data=tmp[1,], aes(x = 0, y = 0.75), size=2, 
#                      label = lm_eqn(x=tmp$m_t_2, y=tmp$m_t, 
#                                     df=tmp), parse = TRUE, vjust=0, hjust=0) 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_CAMDAC_decon_capped_tn2_DMP.pdf"), 
       plot = p1, device = "pdf", units = "in", width = 3.4, height = 2.8)

tmp <- df_phased_methylation[abs(m_diff_tn) > 0.2]
p1 <- ggplot(tmp, aes(x=m_t_2_corr, y=m_t_corr)) + theme_classic() + 
  geom_point(aes(col=vaf_diff), size=1, shape=16, alpha=0.25) + #, col="navyblue") +
  ylim(0,1) + xlim(0,1) + scale_color_gradient(name="VAF\ndiff", low = "blue2", high = "red2") + 
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
            label=paste0("n = ", nrow(tmp), "\ncor = ", 
                         round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]),  digits = 2)), hjust=0)
# p1 <- p1 + geom_text(data=tmp[1,], aes(x = 0, y = 0.75), size=2, 
#                      label = lm_eqn(x=tmp$m_t_2, y=tmp$m_t, 
#                                     df=tmp), parse = TRUE, vjust=0, hjust=0) 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_CAMDAC_decon_capped_tn_DMP.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

p1 <- ggplot(tmp, aes(x=m_t_2_corr-m_t_corr, y=..density..)) + 
  geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() + 
  geom_text(data=df_phased_methylation[1,], aes(x = -max_minus+0.5, y = 3), size=2, 
            label=paste0("n = ", nrow(tmp),  "\ncor = ", 
                         round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]), digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_CAMDAC_decon_density_capped.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

df <- data.frame(y = c(tmp[,(m_t_2_corr)], tmp[,(m_t_corr)]),
                 x = factor(c(rep(0, nrow(tmp)),rep(1, nrow(tmp))), levels = c(0,1), ordered = TRUE), 
                 flag = c(rep("SNV", nrow(tmp)),rep("CAMDAC", nrow(tmp))))
my_comparisons <- list(factor(c(0, 1), levels = c(0,1), ordered = TRUE))

df_phased_methylation[,median(abs(m_t_2-m_t), na.rm=T)] 
mads <- df_phased_methylation[,mad(SNV_VAF_methylation-SNV_VAF_WGS, constant = 3, na.rm=T)]
tmp <- df_phased_methylation[abs(SNV_VAF_methylation-SNV_VAF_WGS) < mads,]

p1 <- ggplot(tmp, aes(x=m_t_2-m_t, y=..density..)) + 
  geom_density(alpha=0.25, col="navyblue", fill="navyblue") + theme_classic() + 
  geom_text(data=df_phased_methylation[1,], aes(x = -1, y = 3), size=2, 
            label=paste0("n = ", nrow(tmp),  "\ncor = ", 
                         round2(cor(tmp[, m_t], tmp[, m_t_2]), digits = 2)), hjust=0)
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/",
                          "SNV_decon_vs_CAMDAC_decon_density_filtered.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

p1 <- ggplot(tmp, aes(x=m_t_2, y=m_t)) + theme_classic() +
  geom_point(size=0.5, shape=16, alpha=0.25, col="navyblue") + 
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
  geom_text(data=tmp[1,], aes(x = 0, y = 1.12), size=2, 
            label=paste0("n = ", nrow(tmp),  "\ncor = ", 
                         round2(cor(tmp[, m_t], tmp[, m_t_2]), digits = 2)), hjust=0)
# p1 <- p1 + geom_text(data=tmp[1,], aes(x = 0, y = 1.2), size=2, 
#                      label = lm_eqn(x=tmp$m_t_2, y=tmp$m_t, df=tmp), parse = TRUE, vjust=0, hjust=0) 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_CAMDAC_filtered.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3, height =2.8)

p1 <- ggplot(tmp, aes(x=m_t_2_corr, y=m_t_corr)) + theme_classic() + 
  geom_point(size=1, shape=16, alpha=0.25, col="black", fill="black") + 
  ylim(0,1) + xlim(0,1) + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_abline(slope=1, intercept = 0, size=0.5, linetype="dashed", col="navyblue") +
  geom_text(data=tmp[1,], aes(x = 0, y = 0.85), size=2, 
            label=paste0("n = ", nrow(tmp), "\ncor = ", 
                         round2(cor(tmp[, m_t_corr], tmp[, m_t_2_corr]),  digits = 2)), hjust=0)
# p1 <- p1 + geom_text(data=tmp[1,], aes(x = 0, y = 0.75), size=2, 
#                      label = lm_eqn(x=tmp$m_t_2, y=tmp$m_t, 
#                                     df=tmp), parse = TRUE, vjust=0, hjust=0) 
ggsave(filename =  paste0("~/rdeocd/txMeth_analysis/SNV_decon/SNV_decon_vs_CAMDAC_decon_capped_filtered.pdf"), 
       plot=p1, device = "pdf", units = "in", width = 3.4, height =2.8)

rm(list=ls(pattern = "^p[0-9]|^df$|tmp|mads|df_phased_methylation|minus|plus"))


