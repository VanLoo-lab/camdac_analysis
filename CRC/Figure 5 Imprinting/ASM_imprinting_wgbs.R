### Most updated WGBS imprinting analysis code adapted from Lili's code and using Lili's input files 
### Author: Lesky

setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting")

# Load library
x <- c("GenomicRanges","data.table","stringr","readr","ggplot2","gridExtra")
# library(pryr) # only needed if mem_used() is used to determine object RAM.
invisible(lapply(x, function(y) suppressWarnings(
  suppressMessages(library(y, character.only = TRUE,quietly=TRUE,
                           warn.conflicts = FALSE)))))
rm(x)
library(fst)
library(tidyverse)

# read imprinted genes
df_imprinted_genes <- as.data.table(read_csv("Imprinted_lili_hg38.csv"))
genes <- df_imprinted_genes$Gene


# load cancer genes
# COSMIC
dt_cancer <- data.table(read.table(
  "../cancer_gene_census.csv", 
  header = TRUE, sep = ",", stringsAsFactors = FALSE))
dt_cancer <- dt_cancer[grepl("(CRC)|(colorectal)",Tumour.Types.Somatic.) & !grepl("CCRCC", Tumour.Types.Somatic.),]

# PANCAN
dt_cancer3 <- data.table(read.table(
  "../20200401_pancan_drivers.csv", 
  header = TRUE, sep = ",", stringsAsFactors = FALSE))
cancer_genes = sort(unique(c(dt_cancer$Gene.Symbol, dt_cancer3$gene)))
df_imprinted_genes[, cancer_driver := ifelse(Gene %in% cancer_genes, TRUE, FALSE)]
rm(list=ls(pattern="dt_cancer"))

# set genome build var
build="hg38"


# annotate each CpG to a bin
all_regions_anno <- read_fst("/rsrch4/home/genetics/slai4/R/x86_64-pc-linux-gnu-library/4.1/CAMDAC/pipeline_files/hg38_annotations/hg38_all_regions_annotations.fst", as.data.table=TRUE)
all_regions_anno <- all_regions_anno[!grepl("X|Y", chrom)]
all_regions_anno[, chrom := as.numeric(chrom)]
setkeyv(all_regions_anno, cols=c("chrom", "start", "end"))
setkeyv(df_imprinted_genes, cols=c("chrom", "start", "end"))
all_regions_anno <- 
  all_regions_anno[, .SD, .SDcols=colnames(all_regions_anno)[!grepl("(_ids)|(repeat_)", colnames(all_regions_anno))]]
ov <- foverlaps(all_regions_anno, df_imprinted_genes)
ov <- ov[!is.na(Gene)&!is.na(cluster_id)]

# re-format ov columns
colnames(ov)[grepl("(^start)|(^end)", colnames(ov))] <- c("bin_start_lit", "bin_end_lit")
colnames(ov) <- gsub("i\\.", "", colnames(ov))
colnames(ov)[1] <- "chrom"
ov[, chrom := as.factor(chrom)]

# collapse duplicate ids
ov <- ov[,.(Gene=paste(sort(unique(Gene)), collapse = ","),
            Status=paste(sort(unique(Status)), collapse = ","),
            driver=ifelse(sum(as.numeric(cancer_driver))>0, TRUE, FALSE)), 
         by=.(cluster_id, chrom, start, end, intragenic, CGI, enhancer, promoter, 
              exon, intron, FUTR, TUTR, repeats)]
setkeyv(ov, cols=c("chrom", "start", "end"))



# load patient_info
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
asm_paths <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/N/", PT_ID, ".N.asm_ss_dmp.csv.gz"))
bulk_paths <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/N/", PT_ID, ".N.m.csv.gz")) 


# get consectuvive DMP helper function
fct <- function(x,sgn,min) {
  flag <- c(FALSE, diff(sign(sgn)) > 0)
  y <- rle(x==TRUE & flag==FALSE)
  z<-as.numeric(suppressWarnings(
    max(y$lengths[y$lengths >= min & y$values == TRUE]))) 
  if(!is.finite(z)){z <- as.numeric(NA)}
  return(z)
} 


# ### testing only: check how DMP are called and why ASM_dmp_b are all NA ---
# itersplit <- 1e5
# effect_size <- 0.5
# prob <- 0.99
# dt <- data.table::fread(asm_paths[1])  
# asm_b_diff <- dt[["alt_m"]] - dt[["ref_m"]]
# M_alt <- round(dt[["alt_total_counts_m"]] * dt[["alt_m"]], 0)
# UM_alt <- dt[["alt_total_counts_m"]] - M_alt
# M_ref <- round(dt[["ref_total_counts_m"]] * dt[["ref_m"]], 0)
# UM_ref <- dt[["ref_total_counts_m"]] - M_ref

# # asm_b_diff <- asm_b_diff[1:20]
# # M_alt <- M_alt[1:20]
# # UM_alt <- UM_alt[1:20]
# # M_ref <- M_ref[1:20]
# # UM_ref <- UM_ref[1:20]
# asm_t_diff <- NULL

# # dmp calling pipeline 
# # dmp_call <- CAMDAC:::dmp_call_pipe(asm_b_diff, M_ref, UM_ref, M_alt, UM_alt, asm_t_diff, effect_size, prob, itersplit)

# # step by step
# phypo <- CAMDAC:::calc_prob_dmp(M_ref, UM_ref, M_alt, UM_alt, ncores = ncores, itersplit = itersplit)
# prob_DMP <- data.table::fifelse(asm_b_diff > 0, 1 - phypo, phypo)
# asm_DMP_b <- CAMDAC:::prob_to_call(prob_DMP, asm_b_diff, effect_size, prob)
# dmp_call <- data.table(M_ref, UM_ref, M_alt, UM_alt, prob_DMP, asm_b_diff, asm_DMP_b)
# dt[,12:14] <- NULL
# dt <- cbind(dt, dmp_call)

# prob_to_call <- function(p, mdiff, effect_size = 0.2, prob = 0.99) {
#   data.table::fcase(
#     p >= prob & mdiff >= effect_size, "hyper",
#     p >= prob & mdiff <= (-effect_size), "hypo"
#   )
# }

# # compare with HDI overlap method: much more DMPs with Nana's method
# dt[, is_DM_n12 := (ref_m_hi<alt_m_lo|alt_m_hi<ref_m_lo) & abs(asm_b_diff) > 0.5]
# table(dt$asm_DMP_b, useNA="ifany")
# table(dt$is_DM_n12, useNA="ifany")
# table(dt$asm_DMP_b, dt$is_DM_n12, useNA="ifany")

# #####


# get imprinted bins
dt_IM_regions_list <- list()
for(i in 1:6){
  print(i)
  df_phased_methylation <- data.table::fread(asm_paths[i]) 
  bulk_n <- data.table::fread(bulk_paths[i])
  bulk_n <- bulk_n[, .SD, .SDcols=c("chrom", "start", "end", "m", "m_x_low", "m_x_high")]
  setnames(bulk_n, old = c("m", "m_x_low", "m_x_high"), new = c("m_n", "m_n_low", "m_n_high"))
  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  setkeyv(bulk_n, cols=c("chrom", "start", "end"))
  setkeyv(ov, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, bulk_n, nomatch = 0)
  # identical(df_phased_methylation$start, df_phased_methylation$i.start)
  df_phased_methylation[, `:=`(i.start = NULL, i.end = NULL)] # make sure i.start i.end are identical to start end !!!
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, ov)
  df_phased_methylation <- df_phased_methylation[!is.na(chrom)&!is.na(Gene)]
  setnames(df_phased_methylation, 
           old = c("asm_b_diff", "asm_DMP_b", "ref_m_lo", "ref_m_hi", "alt_m_lo", "alt_m_hi"), 
           new = c("m_diff_n12", "is_DM_n12", "m_n_1_low", "m_n_1_high", "m_n_2_low", "m_n_2_high"))
  # df_phased_methylation$is_DM_n12 <- ifelse(abs(df_phased_methylation$m_diff_n12)>0.5, TRUE, FALSE)
  df_phased_methylation[, is_DM_n12 := (m_n_1_high<m_n_2_low|m_n_2_high<m_n_1_low) & abs(m_diff_n12) > 0.5]
  # annotate imprinted gene loci
  # rename the CpG start and end position
  colnames(df_phased_methylation) <- gsub("i\\.", "CpG_", colnames(df_phased_methylation))
  
  # filter out CpGs with insufficient counts and/or non-CG SNPs
  df_phased_methylation <- df_phased_methylation[!is.na(m_n)]
  df_phased_methylation <- df_phased_methylation[!is.na(m_diff_n12)]
  
  # filter out CpGs with intermetdiate methylation
  df_phased_methylation[, IM := m_n_low>0.2&m_n_high<0.8]
    
  # check imprinted genes
  # df_phased_methylation[, cluster_id := paste(cluster_id, Gene, sep="_")]
  tmp <- df_phased_methylation[, .SD, .SDcols=c("m_diff_n12", "IM", "cluster_id")]
  tmp <- tmp[,.(consec_IM = fct(x = IM, sgn=m_diff_n12, min = 2)),
               by=.(cluster_id)] 
  df_phased_methylation$consec_IM_m_n <- tmp[match(df_phased_methylation$cluster_id, cluster_id), 
                                                .SD, .SDcols=c("consec_IM")] 
  rm(tmp)
    
  # Annotate imprinted genes according to consecutive bulk IM CpGs
  df_phased_methylation[, IM_average := ifelse(is.na(consec_IM_m_n), "not_IM",
                                        ifelse(consec_IM_m_n >= 2, "IM","unknown"))]
    
  # check imprinted genes
  tmp <- df_phased_methylation[, .SD, .SDcols=c("is_DM_n12", "m_diff_n12", "cluster_id")]
  tmp <- tmp[, .(consec_IM = fct(x = (!is.na(is_DM_n12) & is_DM_n12 == TRUE & abs(m_diff_n12) > 0.7), 
                                 sgn = m_diff_n12, min = 2)), by=.(cluster_id)] 
  df_phased_methylation$consec_IM_m_n_phased <- tmp[match(df_phased_methylation$cluster_id, cluster_id), 
                                                  .SD, .SDcols=c("consec_IM")] 
  rm(tmp)
    
  # Annotate imprinted genes according to consecutive phased IM CpGs
  df_phased_methylation[, IM_phased := ifelse(is.na(consec_IM_m_n_phased), "not_IM",
                                              ifelse(consec_IM_m_n_phased >= 2, "IM","unknown"))]
  df_phased_methylation[, flag := !is.na(is_DM_n12) & is_DM_n12 == TRUE & abs(m_diff_n12) > 0.6]
  
  if(sum(!is.na(df_phased_methylation$IM_phased)&df_phased_methylation$IM_phased=="IM")>0){
    # save imprinted gene coordinates
    dt_IM_regions_list[[length(dt_IM_regions_list)+1]] <- 
        df_phased_methylation[!is.na(IM_phased)&IM_phased=="IM", 
                              .(chrom=as.character(chrom[1]),  start=min(CpG_start[flag==TRUE], na.rm=TRUE), end=max(CpG_end[flag==TRUE], na.rm=TRUE), 
                                CpG_counts=ifelse(sum(flag==TRUE)==0, as.numeric(NA), max(which(flag==TRUE)) - min(which(flag==TRUE)) + 1), 
                                consec_IM_m_n=mean(consec_IM_m_n, na.rm=TRUE), consec_IM_m_n_phased=mean(consec_IM_m_n_phased, na.rm=TRUE), 
                                IM_m_n_average = IM_average[1], IM_m_n_phased = IM_phased[1]), 
                            by=.(cluster_id, Gene, Status, driver)]
  }
} ; rm(df_phased_methylation)
dt_IM_regions <- rbindlist(dt_IM_regions_list)
table(dt_IM_regions$IM_m_n_phased) # 50 bins with imprinting; there are overlaps
rm(dt_IM_regions_list)

# extract CAMDAC cluster id
# dt_IM_regions[, id := tstrsplit(cluster_id, "_", keep=1)]
dt_IM_regions[, id := cluster_id]

# regions list
dt_IM_regions <- dt_IM_regions[, .(Gene=paste(sort(unique(Gene)), collapse = ","),
                                   Status=paste(sort(unique(Status)), collapse = ","),
                                   driver=ifelse(sum(as.numeric(driver))>0, TRUE, FALSE),
                                   chromosome=as.character(chrom[1]), start=min(start), end=max(end),
                                   CpGs_bin = ifelse(sum(!is.na(CpG_counts))==0, as.numeric(NA), max(CpG_counts, na.rm=TRUE)),
                                   consec_IM_m_n=ifelse(sum(!is.na(consec_IM_m_n))>0, 
                                                        max(consec_IM_m_n, na.rm=TRUE), as.numeric(NA)),
                                   consec_IM_m_n_phased=ifelse(sum(!is.na(consec_IM_m_n_phased))>0, 
                                                               max(consec_IM_m_n_phased, na.rm=TRUE), as.numeric(NA)),
                                   IM_m_n_average = sum(IM_m_n_average=="IM"&!is.na(IM_m_n_average))>0, 
                                   IM_m_n_phased = sum(IM_m_n_phased=="IM"&!is.na(IM_m_n_phased))>0), 
                  by=.(id)] 
dt_IM_regions <- dt_IM_regions[IM_m_n_phased==TRUE][order(id, Gene),] 

# reformat gene loci 
dt_IM_regions[grepl("GNAS", Gene), Gene := "GNAS,GNASAS"]
dt_IM_regions[grepl("GNAS", Gene), cancer_driver := TRUE]

# prep for overlaps
setkeyv(dt_IM_regions, cols=c("chromosome", "start", "end"))
all_regions_anno[, chrom := as.character(chrom)]
setkeyv(all_regions_anno, cols=c("chrom", "start", "end"))
# annotate each bin
dt_IM_regions_anno <- foverlaps(all_regions_anno[, .SD, .SDcols=colnames(all_regions_anno)[colnames(all_regions_anno)!="cluster_id"]],
                dt_IM_regions[, .SD, .SDcols = c(1, which(colnames(dt_IM_regions)=="Gene"):which(colnames(dt_IM_regions)=="CpGs_bin"))])
dt_IM_regions_anno <- dt_IM_regions_anno[!is.na(Gene)]
rm(all_regions_anno, dt_IM_regions)

# re-format dt_IM_regions_anno
colnames(dt_IM_regions_anno) <- gsub("i\\.", "bin_", colnames(dt_IM_regions_anno))
colnames(dt_IM_regions_anno)[1] <- "chrom"
setkeyv(dt_IM_regions_anno, cols=c("chrom", "start", "end")) # 12 unique genes with 24 bins 

# get normal methylation at imprinted bins
dt_IM_regions_normal_list <- list()

for(i in 1:6){
  print(i)
  df_phased_methylation <- data.table::fread(asm_paths[i]) 
  bulk_n <- data.table::fread(bulk_paths[i])
  bulk_n <- bulk_n[, .SD, .SDcols=c("chrom", "start", "end", "m", "m_x_low", "m_x_high")]
  setnames(bulk_n, old = c("m", "m_x_low", "m_x_high"), new = c("m_n", "m_n_low", "m_n_high"))
  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  setkeyv(bulk_n, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, bulk_n, nomatch = 0)
  df_phased_methylation[, `:=`(i.start = NULL, i.end = NULL)]
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, dt_IM_regions_anno)
  df_phased_methylation <- df_phased_methylation[!is.na(chrom)&!is.na(Gene)]
  setnames(df_phased_methylation, 
           old = c("asm_b_diff", "asm_DMP_b", "ref_m_lo", "ref_m_hi", "alt_m_lo", "alt_m_hi"), 
           new = c("m_diff_n12", "is_DM_n12", "m_n_1_low", "m_n_1_high", "m_n_2_low", "m_n_2_high"))
  
  # annotate imprinted gene loci
  colnames(df_phased_methylation) <- gsub("i\\.", "CpG_", colnames(df_phased_methylation))
  
  # filter out CpGs with insufficient counts and/or non-CG SNPs
  df_phased_methylation <- df_phased_methylation[!is.na(m_n)]
  df_phased_methylation <- df_phased_methylation[!is.na(m_diff_n12)]
  
  # only keep loci with 3 phased CpGs or more
  df_phased_methylation[, CpGs := length(chrom), by=.(id)]
  df_phased_methylation <- df_phased_methylation[CpGs >= 3,]

  # filter out CpGs with intermetdiate methylation
  df_phased_methylation[, IM := m_n_low>0.2&m_n_high<0.8]

  # correct DM threshold
  df_phased_methylation[, is_DM_n12 := (m_n_1_high<m_n_2_low|m_n_2_high<m_n_1_low) & abs(m_diff_n12) > 0.5]
  
  # check imprinted genes
  tmp <- df_phased_methylation[, .SD, .SDcols=c("m_diff_n12", "IM", "id")]
  tmp <- tmp[,.(consec_IM = fct(x = IM, sgn=m_diff_n12, min = 2)),
             by=.(id)] 
  df_phased_methylation$consec_IM_m_n <- 
    tmp[match(df_phased_methylation$id, id),.SD, .SDcols=c("consec_IM")] 
  rm(tmp)
  
  # Annotate imprinted genes according to consecutive bulk IM CpGs
  df_phased_methylation[, IM_average := ifelse(is.na(consec_IM_m_n), "NA",
                                        ifelse(consec_IM_m_n >= 2, "IM","not_IM"))]
  
  # check imprinted genes
  tmp <- df_phased_methylation[, .SD, .SDcols=c("is_DM_n12", "m_diff_n12", "id")]
  tmp <- tmp[, .(consec_IM = fct(x = (!is.na(is_DM_n12) & is_DM_n12 == TRUE), 
                                 sgn = m_diff_n12, min = 2)), by=.(id)] 
  df_phased_methylation$consec_IM_m_n_phased <- 
    tmp[match(df_phased_methylation$id, id), .SD, .SDcols=c("consec_IM")] 
  rm(tmp)
  
  # Annotate imprinted genes according to consecutive phased IM CpGs
  df_phased_methylation[, IM_phased := ifelse(is.na(consec_IM_m_n_phased), "NA",
                                              ifelse(consec_IM_m_n_phased >= 2, "IM","not_IM"))]
  
  # save imprinted gene coordinates
  df_phased_methylation[, patient := PT_ID[i]]
  dt_IM_regions_normal_list[[length(dt_IM_regions_normal_list)+1]] <- 
    df_phased_methylation[, .(chrom=as.character(chrom[1]), start=start[1], end=end[1],
                              bin_start=bin_start[1], bin_end=bin_end[1],
                              consec_IM_m_n=mean(consec_IM_m_n, na.rm=TRUE),
                              consec_IM_m_n_phased=mean(consec_IM_m_n_phased, na.rm=TRUE), 
                              IM_m_n_average = IM_average[1], IM_m_n_phased = IM_phased[1]), 
                          by=.(patient, id, Gene, Status, driver, 
                               CpGs_bin, CpGs,
                               intragenic, CGI, enhancer, promoter, exon, 
                               intron, FUTR, TUTR, repeats)]

  # Load bulk methylation information
  df_bulk <- data.table::fread(bulk_paths[i])
  df_bulk <- df_bulk[, .(chrom, start, end, m, m_x_low, m_x_high)]
  setnames(df_bulk, old = c("m", "m_x_low", "m_x_high"), new = c("m_n", "m_n_low", "m_n_high"))  
  
  for(ID in unique(df_phased_methylation$id)){
    # extact gene
    tmp <- df_phased_methylation[id==ID,]
    gene <- unique(tmp$Gene)

    # extact coordinates
    coords <- dt_IM_regions_anno[Gene==gene&id==ID, .(chrom, "start"=bin_start, "end"=bin_end)] 
    coords$anno<-coords$start+(1/6)*(coords$end-coords$start)
    
    # re-format gene
    gene <- gsub(",", "_", gene)
    print(paste0(gene, "_", ID))
    
    # Get methylation and stats info or each CpG overlapping with the promoter regions
    df_anno <- df_bulk[chrom==as.character(coords$chrom[1])&
                       start >= min(c(coords$start - 400, tmp$start)) & 
                       end <= max(c(coords$end + 400, tmp$end)),]
      
    p <- ggplot(tmp) +  theme_classic() + 
      ggtitle(paste(gene, "_", ID, " (", unique(tmp$Status),")", sep="")) +
      scale_y_continuous(name = "methylation rates", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), limits=c(-0.15,1.15))+
      coord_cartesian(xlim=c(min(df_anno$start)-100, max(df_anno$end)+100), expand = FALSE) +
      xlab("Chromosome Position (bp)") +
      theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) + 
      scale_color_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))+
      scale_fill_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))
    
    p1 <- p + geom_point(data=df_anno, aes(x=start,y=m_n),col="grey10")+
      geom_ribbon(data=df_anno, aes(x=start,ymin=m_n_low,ymax=m_n_high), fill="grey10", alpha=.2)+
      geom_point(data=tmp, aes(x=CpG_start,y=ref_m,col="ref"))+
      geom_point(data=tmp, aes(x=CpG_start,y=alt_m,col="alt"))
    
    plot_outdir <- paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_WGBS/", PT_ID[i], "/N/")
    if (!dir.exists(plot_outdir)) {
      dir.create(plot_outdir, recursive = TRUE)
      }
    
    ggsave(filename =  paste0(plot_outdir, gene, "_", ID,"_plot.pdf"), 
           plot=p1, device = "pdf", units = "in", width = 5, height = 3)
    rm(p1, p, df_anno, tmp, coords, p1)
  }
  rm(df_bulk, df_phased_methylation)
} 
dt_IM_normals <- rbindlist(dt_IM_regions_normal_list)
save(dt_IM_normals, dt_IM_regions_anno, file="ASM/wgbs_normal_imprinting_20240704.RData")

load("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM/wgbs_normal_imprinting_20240704.RData")
# number of imprinted loci
table(dt_IM_normals$IM_m_n_phased)
df_gene_patient <- as.data.frame(table(dt_IM_normals$patient, dt_IM_normals$Gene))
df_gene_patient$Freq <- ifelse(df_gene_patient$Freq > 0, 1, 0)
library(dplyr)
df_gene_patient %>%
  group_by(Var1) %>%
  summarise(Freq = sum(Freq))


# Imprinting in tumors -------
dt_IM_tumours_list <- list()
asm_paths_t <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/T/", PT_ID, ".T.asm_meth_pure.csv.gz"))
bulk_paths_t <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/T/", PT_ID, ".T.pure.csv.gz"))

for(i in 1:6){
  print(i)
  df_phased_methylation <- data.table::fread(asm_paths_t[i]) 
  setnames(df_phased_methylation, 
           old = c("ref_m_t", "alt_m_t", "ref_m_t_low", "ref_m_t_high", "alt_m_t_low", "alt_m_t_high", "ref_m_i", "alt_m_i", "ref_CN", "alt_CN"), 
           new = c("m_t_1", "m_t_2", "m_t_1_low", "m_t_1_high", "m_t_2_low", "m_t_2_high", "m_n_1", "m_n_2", "n1", "n2")) 

  # Lesky's comment: add bulk, normal, asm normal methylation rate 
  bulk_t <- data.table::fread(bulk_paths_t[i])
  bulk_t <- bulk_t[, .SD, .SDcols=c("chrom", "start", "end", "m_t", "m_t_low", "m_t_high", "m_n")]
  setnames(bulk_t, old = "m_t", new = "m_t_corr") 
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  setkeyv(bulk_t, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, bulk_t, nomatch = 0)
  df_phased_methylation[, `:=`(i.start = NULL, i.end = NULL)]

  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, dt_IM_regions_anno)
  df_phased_methylation <- df_phased_methylation[(!is.na(m_t_1)|!is.na(m_t_2))&!is.na(Gene)]
  
  # filter out CpGs with insufficient counts and/or non-CG SNPs
  # Lesky's comment: check this label in RRBS 
  # df_phased_methylation <- df_phased_methylation[is.na(label)]
  
  # remove duplicate colnames
  colnames(df_phased_methylation) <- gsub("i\\.", "CpG_", colnames(df_phased_methylation))
  df_phased_methylation <- df_phased_methylation[!is.na(Gene)][order(id, Gene, CpG_start, CpG_end)]
  
  # filter out bins IM and validated by SNP in normal
  tmp <- dt_IM_normals[patient==PT_ID[i]&IM_m_n_phased=="IM"]
  df_phased_methylation <- df_phased_methylation[id %in% tmp[CpGs>2,unique(id)]]
  df_phased_methylation[, as.character(c("CpGs_n", "IM_m_n_average" , "consec_IM_m_n", 
                                          "IM_m_n_phased", "consec_IM_m_n_phased")) := 
                          as.list(tmp[match(df_phased_methylation$id, id), 
                                      .(CpGs, IM_m_n_average , consec_IM_m_n, 
                                        IM_m_n_phased, consec_IM_m_n_phased)])]
  
  
  # filter out CpGs with intermetdiate methylation
  df_phased_methylation[, IM_t_average := ifelse(is.na(m_t_1)&is.na(m_t_2), "no_info", 
                                          ifelse(m_t_low>0.2&m_t_high<0.8, "IM", "not_IM"))]

  # check imprinted genes
  df_phased_methylation[, m_diff_t12 := abs(m_t_1 - m_t_2)] 
  tmp <- df_phased_methylation[, .SD, .SDcols=c("m_diff_t12", "IM_t_average", "id")]
  tmp <- tmp[,.(consec_IM = fct(x = IM_t_average=="IM", sgn=m_diff_t12, min = 2)), by=.(id)]
  df_phased_methylation$consec_IM_m_t_average <-
      tmp[match(df_phased_methylation$id, id), .SD, .SDcols=c("consec_IM")]
  rm(tmp)

  # Annotate imprinted genes according to consecutive bulk IM CpGs
  df_phased_methylation[, IM_m_t_average := ifelse(is.na(m_t_1)&is.na(m_t_2), "no_info", 
                                            ifelse(!is.na(consec_IM_m_t_average) &
                                                    consec_IM_m_t_average >= 2, "IM","not_IM"))]
  # df_phased_methylation[, is_DM_t12 := abs(m_diff_t12) > 0.5] 
  # Lesky's comment: check HDI 
  df_phased_methylation[, is_DM_t12 := (m_t_1_high < m_t_2_low | m_t_2_high < m_t_1_low) & 
                          abs(m_diff_t12) > 0.5]
  
  # check imprinted genes
  tmp <- df_phased_methylation[, .SD, .SDcols=c("is_DM_t12", "m_diff_t12", "id")]
  tmp <- tmp[, .(consec_IM = fct(x = (!is.na(is_DM_t12) & is_DM_t12 == TRUE & abs(m_diff_t12) > 0.5), 
                                  sgn = m_diff_t12, min = 2)), by=.(id)] 
  df_phased_methylation$consec_IM_m_t_phased <- 
      tmp[match(df_phased_methylation$id, id), .SD, .SDcols=c("consec_IM")] 
  rm(tmp)
  
  # Annotate imprinted genes according to consecutive DM alleles
  df_phased_methylation[, IM_m_t_phased := ifelse(is.na(m_t_1) & is.na(m_t_2), "no_info", 
                                            ifelse(!is.na(consec_IM_m_t_phased) & 
                                                    consec_IM_m_t_phased >= 2, "IM",
                                                    "not_IM"))]
    
  # save imprinted gene coordinates
  df_phased_methylation[, patient := PT_ID[i]]
  df_phased_methylation[, LOH := ifelse((n1==0|n2==0), 1, 0)]
  df_phased_methylation[, prob_decon := ifelse(m_t_corr > m_n & n1 == 0, m_t_2 + 0.1 < m_n,
                                        ifelse(m_t_corr < m_n & n1 == 0, m_t_2 - 0.1 > m_n,
                                        ifelse(m_t_corr > m_n & n2 == 0, m_t_1 + 0.1 < m_n,
                                        ifelse(m_t_corr < m_n & n2 == 0, m_t_1 - 0.1 > m_n, as.logical(NA)))))]

  # summarise info by gene
  dt_info <- df_phased_methylation[, .(chrom=as.character(chrom[1]), 
                                         start=min(start), end=max(end), 
                                         CpGs_t=length(chrom),
                                         LOH=sum(LOH, na.rm=TRUE)>=1, 
                                         n1 = mean(n1, na.rm=TRUE), 
                                         n2 = mean(n2, na.rm=TRUE),
                                         m_diff_n12=abs(mean(m_n_1 - m_n_2, na.rm=TRUE)), 
                                         m_diff_t12=abs(mean(m_t_1 - m_t_2, na.rm=TRUE)),
                                         m_n_1=mean(m_n_1, na.rm=TRUE),
                                         m_n_2=mean(m_n_2, na.rm=TRUE),
                                         m_t_1=mean(m_t_1, na.rm=TRUE),
                                         m_t_2=mean(m_t_2, na.rm=TRUE),
                                         m_t=mean(m_t_corr, na.rm=TRUE),
                                         consec_IM_m_n=mean(consec_IM_m_n, na.rm=TRUE),
                                         consec_IM_m_n_phased=mean(consec_IM_m_n_phased, na.rm=TRUE),
                                         consec_IM_m_t_average=mean(consec_IM_m_t_average, na.rm=TRUE),
                                         consec_IM_m_t_phased=mean(consec_IM_m_t_phased, na.rm=TRUE),
                                         prob_decon=sum(prob_decon, na.rm=TRUE)>0), 
                                       by=.(patient, id, 
                                            Gene, Status, 
                                            CpGs_bin, driver, CpGs_n, 
                                            IM_m_n_average, IM_m_n_phased,
                                            IM_m_t_average, IM_m_t_phased)]

    dt_info[CpGs_t < 3, IM_m_t_phased := "no_info"]
    dt_info[CpGs_t < 3, IM_m_t_average := "no_info"]
    dt_info[, flag := ifelse(IM_m_t_phased == "no_info", "no_info", 
                      ifelse(LOH == TRUE, "LOH", 
                      ifelse(IM_m_t_phased == "IM", "Imprinted_phased",
                      ifelse(m_diff_t12 >= 0.5, "Imprinted_maybe", 
                      ifelse(m_diff_t12 <= 0.15 & IM_m_t_phased == "not_IM", "LOI", 
                      ifelse(m_diff_t12 < 0.5, "partial_LOI", "investigate"))))))]
  
    dt_IM_tumours_list[[length(dt_IM_tumours_list)+1]] <- dt_info

    # plot imprinted genes in tumor after CAMDAC ----------------
    bulk_t <- data.table::fread(bulk_paths_t[i])
    df_bulk <- bulk_t[, .SD, .SDcols=c("chrom", "start", "end",
                                                 "m_n", "cov_n", "m_x_low_n", "m_x_high_n",
                                                 "m", "cov", "m_x_low", "m_x_high",
                                                 "m_t", "cov_t", "m_t_low", "m_t_high")]
    setnames(df_bulk, old = c("m_x_low_n", "m_x_high_n", "m", "cov", "m_x_low", "m_x_high", "m_t"), new = c("m_n_low", "m_n_high", "m_b", "cov_b", "m_b_low", "m_b_high", "m_t_corr")) 
    for(ID in unique(df_phased_methylation$id)){
      print(ID)
      # extact gene
      tmp <- df_phased_methylation[id==ID,]
      gene <- unique(tmp$Gene)

      # extact coordinates
      coords <- dt_IM_regions_anno[id==ID, .SD,
                      .SDcols=c("chrom", "start", "end")]
      coords$anno<-coords$start+(1/6)*(coords$end-coords$start)
      
      # Get methylation and stats info or each CpG overlapping with the promoter regions
      df_anno <- df_bulk[chrom==as.character(coords$chrom[1])&
                            start >= min(c(coords$start-400, tmp$start)) &
                            end <= max(c(coords$end+ 400, tmp$end)),]
      df_anno[, m_t_low := ifelse(m_t_low <0 , 0, ifelse(m_t_low>1, 1, m_t_low))]
      df_anno[, m_t_high := ifelse(m_t_high<0 , 0, ifelse(m_t_high>1, 1, m_t_high))]
      
      p <- ggplot(tmp) +  theme_classic() +
        ggtitle(paste(gene, "_", ID, " (", unique(tmp$Status), ") ASCN= ",
                      unique(tmp$n1), " + ", unique(tmp$n2), sep="")) +
        scale_y_continuous(name = "methylation rates", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), limits=c(-0.15,1.15))+
        coord_cartesian(xlim=c(min(df_anno$start)-100, max(df_anno$end)+100), expand = FALSE) +
        xlab("Chromosome Position (bp)") +
        theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) +
        scale_color_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))+
        scale_fill_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))
      
      p1 <- p + geom_point(data=df_anno, aes(x=start,y=m_n),col="grey10")+
          geom_ribbon(data=df_anno, aes(x=start,ymin=m_n_low,ymax=m_n_high), fill="grey10", alpha=.2)+
          geom_point(data=tmp, aes(x=CpG_start,y=m_n_1,col="ref"))+
          geom_point(data=tmp, aes(x=CpG_start,y=m_n_2,col="alt"))
      
      p2 <- p + geom_point(data=df_anno, aes(x=start,y=m_t_corr),col="grey10")+
          geom_ribbon(data=df_anno, aes(x=start,ymin=m_t_low,ymax=m_t_high), fill="grey10", alpha=.2)+
          geom_point(data=tmp, aes(x=CpG_start,y=m_t_1,col="ref"))+
          geom_point(data=tmp, aes(x=CpG_start,y=m_t_2,col="alt"))

      plot_outdir <- paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_WGBS/", PT_ID[i], "/T/")
      if (!dir.exists(plot_outdir)) {
        dir.create(plot_outdir, recursive = TRUE)
      }
    
      my_layout <- rbind(1,2)
      gd <- arrangeGrob(grobs = list(p1,p2), layout_matrix = my_layout)
      ggsave(filename =  paste0(plot_outdir, gene, "_", ID,"_plot.pdf"), 
           plot=gd, device = "pdf", units = "in", width = 6.5, height = 5.5) 

      rm(p1, p2, my_layout, df_anno)
      
    }

    rm(dt_info)
} ; rm(df_phased_methylation)

dt_IM_tumours <- rbindlist(dt_IM_tumours_list)
save(dt_IM_tumours, file = "ASM/wgbs_tumor_imprinting_20240704.RData")
rm(dt_IM_tumours_list)


# check imprinting in normals
load("ASM/wgbs_normal_imprinting_20240704.RData")
table(dt_IM_normals$Gene, dt_IM_normals$patient)

# check LOI 
load("ASM/wgbs_tumor_imprinting_20240704.RData")
table(dt_IM_tumours$flag)
table(dt_IM_tumours$Gene)
dt_IM_tumours$loss_type <- ifelse(dt_IM_tumours$m_t>0.5, "hyper", "hypo")

table(dt_IM_tumours$Gene[dt_IM_tumours$flag %in% c("Imprinted_phased", "Imprinted_maybe")], dt_IM_tumours$patient[dt_IM_tumours$flag %in% c("Imprinted_phased", "Imprinted_maybe")])
table(dt_IM_tumours$Gene[dt_IM_tumours$loss_type=="hyper" & dt_IM_tumours$flag=="LOH"], dt_IM_tumours$patient[dt_IM_tumours$loss_type=="hyper" & dt_IM_tumours$flag=="LOH"])
table(dt_IM_tumours$Gene[dt_IM_tumours$loss_type=="hypo" & dt_IM_tumours$flag=="LOH"], dt_IM_tumours$patient[dt_IM_tumours$loss_type=="hypo" & dt_IM_tumours$flag=="LOH"])
table(dt_IM_tumours$Gene[dt_IM_tumours$loss_type=="hyper" & dt_IM_tumours$flag %in% c("LOI", "partial_LOI")], dt_IM_tumours$patient[dt_IM_tumours$loss_type=="hyper" & dt_IM_tumours$flag %in% c("LOI", "partial_LOI")])
table(dt_IM_tumours$Gene[dt_IM_tumours$loss_type=="hypo" & dt_IM_tumours$flag %in% c("LOI", "partial_LOI")], dt_IM_tumours$patient[dt_IM_tumours$loss_type=="hypo" & dt_IM_tumours$flag %in% c("LOI", "partial_LOI")])

# check number of phased CpGs
library(data.table)
PT_ID <- c("LORC-02", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
# asm_paths <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/*/", PT_ID, "/Methylation/N/", PT_ID, ".N.asm_ss_dmp.csv.gz"))
asm_paths_t <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/06.CAMDAC/", PT_ID, "/Methylation/T/", PT_ID, ".T.asm_meth.csv.gz"))
for (f in asm_paths_t){
  dt <- fread(f)
  print(nrow(dt))
}

asm_paths_t <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/moonshot_wgbs_rerun/06.CAMDAC/", PT_ID, "/Methylation/T/", PT_ID, ".T.asm_meth_pure.csv.gz"))
for (f in asm_paths_t){
  dt <- fread(f)
  print(nrow(dt))
}
