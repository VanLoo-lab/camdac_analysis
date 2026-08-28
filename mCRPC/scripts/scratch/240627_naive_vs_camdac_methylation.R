#' @description Comparing Bismark copy numbers and purity values across 
#' TRACERx
#' 
#' @author Lili Larose Cadieux

# Load library
x <- c("readr", "GenomicRanges", "data.table", "ggplot2", "gridExtra", "fst")
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

################################################################################
# load all normal samples
################################################################################

# set parent directory
tx_dir <- "~/rdeocd/txMeth/"
work_dir <- "~/data/CAMDAC_paper/figure3"
setwd(work_dir)

# create output directory
if(!dir.exists("outputs")) dir.create("outputs")

# load patient_info
f_nm <- paste0(tx_dir,
               "/tx100clinicalinfo_RRBS_published_status",
               "_immune_clusters_purity_ploidy_CpG_cov.txt")
patient_info <- as.data.table(read.delim(f_nm, header = TRUE))
rm(f_nm)

# remove lymph nodes
patient_info <- patient_info[!grepl("LN", Sample),]

# get patient ids
patient_ids <- unique(patient_info$Patient)

# load/save patient-wise methylation information
normal_lists <- list()
for(patient_id in patient_ids){
  print(patient_id)
  
  if (patient_id != "LTX164") {
    path_patient <- paste0(tx_dir, patient_id, "/",sep="")
    f_nm <-
      paste0(path_patient, "/Allelecounts/N/",
             patient_id, ".N.SNPs.CpGs.all.sorted.RData")
    load(f_nm)
    
    # convert to data.table
    gr_normal$mq <- NULL
    dt <- data.table(data.frame(gr_normal, stringAsFactors=FALSE))
    
    # subset to 100k loci 
    set.seed(1234)
    dt_subset <- dt[sample(seq_len(nrow(dt)), 2E5)]
    rm(dt)

    #Get SNP in ref:alt format 
    dt_subset[, SNP := ifelse(!is.na(ref), paste(ref, alt, sep = ""), NA)]
    dt_subset[, offset := ifelse(width == 4, 1, 0)]
    
    # Annotate SNP type
    dt_subset[, type := factor(ifelse(BAF < 0.15,"Homozygous",
                        ifelse(BAF > 0.85,"Homozygous",
                        ifelse(BAF >= 0.15 & BAF <= 0.85, "Heterozygous",
                               NA))), 
                        levels = c("Homozygous", "Heterozygous"))]

    # Extract columns of interest
    cols <- c("width", "start", "end", "offset",
              "SNP", "POS", "BAF", "type",
              "Tr", "Af", "Tf", "Ar",
              "CGf", "CGr", "TGf", "CAr")
    dt_subset_cols <- dt_subset[, .SD, .SDcols=cols]
    
    # subset to CpGs
    dt_subset_cols <- dt_subset_cols[width > 1, ]
    normal_lists[[length(normal_lists)+1]] <- dt_subset_cols
    
    rm(dt, dt_subset, gr_normal)
  }
}
normals <- rbindlist(normal_lists)
save(normals, file = file.path("outputs", "combined_normals.Rdata"))
rm(normal_lists)

# Counts methylated for allele-specific methylation i.e. excluding CG destroying SNPs.
normals[, M_CAMDAC := 
               ifelse(is.na(SNP), CGf + CGr,
               ifelse(start+offset == POS & 
                        grepl(pattern = "([CT][TC])",SNP), CGr,
               ifelse(end-offset == POS & 
                        grepl(pattern = "([GA][AG])",SNP), CGf,
                      CGf + CGr)))]
normals[, M_NAIVE := CGf + CGr] #i.e. Cf and Gr

# TAs are not being counted properly
# count alternate as unmethylated

# Get the total methylation counts --> CG-destroying SNP allele does not contribute. 
# Will be NA if the CG-destroying allele at CG/CCGG is the only one present
normals[, UM_CAMDAC := 
          ifelse(is.na(SNP),
            TGf + CAr,
          ifelse(start + offset == POS &
            grepl(pattern = "([CT][TC])", SNP),
              CAr,
          ifelse(end - offset == POS &
            grepl(pattern = "([GA][AG])", SNP),
              TGf,
              TGf + CAr)))]
normals[, UM_NAIVE :=
          ifelse(is.na(SNP),
            TGf + CAr, #i.e. Cf and Ar
          ifelse(
            (start + offset == POS &
              grepl(pattern = "(CA|CG|AC|GC)", SNP)),
            TGf + CAr + Ar, # Ar = AAr 
          ifelse(
            (end - offset == POS &
              grepl(pattern = "(TG|CG|GT|GC)", SNP)),
            TGf + CAr + Tf, # Tf = TTf 
          ifelse(start + offset == POS &
            grepl(pattern = "(CT|TC)", SNP),
            TGf + CAr + Tr, # Tr representing TAr
          ifelse(end - offset == POS &
            grepl(pattern = "(GA|AG)", SNP),
            TGf + CAr + Af, # Af representing TAf
          NA)))))]

# Get the counts unmethylated
normals[, "total_counts_m_CAMDAC" := UM_CAMDAC + M_CAMDAC]
normals[, "total_counts_m_NAIVE" := UM_NAIVE + M_NAIVE]

# Remove low coverage loci
normals <- normals[total_counts_m_CAMDAC >= 10 & total_counts_m_NAIVE >= 10,]

# Get the counts unmethylated
normals[, m_CAMDAC := M_CAMDAC / total_counts_m_CAMDAC]
normals[, m_NAIVE := M_NAIVE / total_counts_m_NAIVE]
normals[, offset := NULL]

# get correlations between strand specific and total counts at non-polymorphic CpGs
normals[is.na(SNP) & !is.na(CGr / (CGr + CAr)) & !is.na(CGf / (CGf + TGf)), 
  cor(CGr / (CGr + CAr), CGf / (CGf + TGf))] 
#[1] 0.9577482
normals[is.na(SNP) & !is.na(m_NAIVE) & !is.na(m_CAMDAC),
  cor(m_CAMDAC, m_NAIVE)]
#[1] 0.993544
normals[!is.na(SNP) & type == "Heterozygous" & !is.na(m_NAIVE) & !is.na(m_CAMDAC),
  cor(m_CAMDAC, m_NAIVE)]
#[1] 0.9485323
normals[!is.na(SNP) & type == "Heterozygous" & !is.na(m_NAIVE) & !is.na(m_CAMDAC) &
  grepl("(CT)|(TC)|(GA)|(AG)", SNP), 
  cor(m_CAMDAC, m_NAIVE)]
#[1] 0.9307164

# Set color flag for SNP type
normals[, SNP_type := ifelse(is.na(SNP) | type == "Homozygous", "None", 
  ifelse(grepl("(CT)|(TC)|(GA)|(AG)", SNP), "CT", "CAorCG"))]
cols <- c("CT" = "darkblue", "CAorCG" = "gold3", "None" = "grey50")

# create empty plotting square
empty <- ggplot() + geom_point(aes(1, 1), colour = "white")+
         theme(axis.ticks = element_blank(), panel.background = element_blank(), 
               axis.text.x = element_blank(), axis.text.y = element_blank(),           
               axis.title.x = element_blank(), axis.title.y = element_blank())

# Plot the naive vs. CAMDAC methylation rate across a random sample of 3 million CpGs
set.seed(1234)
tmp <- normals[sample(seq_len(nrow(normals)), 3E6)]
d_top <- ggplot(tmp) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type), # col="red", 
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
scatter <- ggplot(tmp) + theme_classic() +
           geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
            alpha = 0.1, shape = 16, size = 1) +
           scale_color_manual(values = cols) +
           theme(legend.position = "none")
d_right <- ggplot(tmp) + theme_classic() + 
           geom_histogram(
            aes(x = m_CAMDAC, y = ..count.., col = SNP_type), # col="blue", 
            fill=NA,  binwidth=0.04, alpha=0.1) +
           scale_color_manual(values = cols) +
           coord_flip() +
           theme(axis.text.x = element_text(angle = 270),
                legend.position = "none")
gd <-
  grid.arrange(d_top, empty, scatter, d_right, 
               ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
ggsave(filename = file.path("outputs", "normal_methylation.png"), 
       plot = gd, device = "png", width=5, height=6)
ggsave(filename = file.path("outputs", "normal_methylation.pdf"), 
       plot = gd, device = "pdf", width=5, height=6)

# subset to all heterozygous SNPs
tmp <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7]
d_top <- ggplot(tmp) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
scatter <- ggplot(tmp) + theme_classic() +
           geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
            alpha = 0.1, shape = 16, size = 1) +
           scale_color_manual(values = cols) +
           theme(legend.position = "none")
d_right <- ggplot(tmp) + theme_classic() + 
           geom_histogram(
            aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
            fill=NA,  binwidth=0.04, alpha=0.1) +
           scale_color_manual(values = cols) +
           coord_flip() +
           theme(axis.text.x = element_text(angle = 270),
                legend.position = "none")
gd <-
  grid.arrange(d_top, empty, scatter, d_right, 
               ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
ggsave(filename = file.path("outputs", "normal_methylation_SNP.png"), 
       plot = gd, device = "png", width=5, height=5)
ggsave(filename = file.path("outputs", "normal_methylation_SNP.pdf"), 
       plot = gd, device = "pdf", width=5, height=5)

# subset to C>T nd G>A SNPs
tmp2 <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7 &
  grepl("(CT)|(TC)|(GA)|(AG)", SNP)]
d_top <- ggplot(tmp2) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
scatter <- ggplot(tmp2) + theme_classic() +
           geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
            alpha = 0.1, shape = 16, size = 1) +
           scale_color_manual(values = cols) +
           theme(legend.position = "none")
d_right <- ggplot(tmp2) + theme_classic() + 
           geom_histogram(
            aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
            fill=NA,  binwidth=0.04, alpha=0.1) +
           scale_color_manual(values = cols) +
           coord_flip() +
           theme(axis.text.x = element_text(angle = 270),
                legend.position = "none")
gd <-
  grid.arrange(d_top, empty, scatter, d_right, 
               ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
ggsave(filename = file.path("outputs", "normal_methylation_SNP_CT.png"), 
       plot = gd, device = "png", width=5, height=5)
ggsave(filename = file.path("outputs", "normal_methylation_SNP_CT.pdf"), 
       plot = gd, device = "pdf", width=5, height=5)

# subset to C>A nd C>G SNPs
tmp2 <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7 &
  grepl("(CA)|(AC)|(GC)|(CG)", SNP)]
d_top <- ggplot(tmp2) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
scatter <- ggplot(tmp2) + theme_classic() +
           geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
            alpha = 0.1, shape = 16, size = 1) +
           scale_color_manual(values = cols) +
           theme(legend.position = "none")
d_right <- ggplot(tmp2) + theme_classic() + 
           geom_histogram(
            aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
            fill=NA,  binwidth=0.04, alpha=0.1) +
           scale_color_manual(values = cols) +
           coord_flip() +
           theme(axis.text.x = element_text(angle = 270),
                legend.position = "none")
gd <-
  grid.arrange(d_top, empty, scatter, d_right, 
               ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
ggsave(filename = file.path("outputs", "normal_methylation_SNP_CAorCG.png"), 
       plot = gd, device = "png", width=5, height=5)
ggsave(filename = file.path("outputs", "normal_methylation_SNP_CAorCG.pdf"), 
       plot = gd, device = "pdf", width=5, height=5)

# clean-up
rm(d_top, empty, d_right, gd, scatter, normals, tmp, tmp2)

################################################################################
# load all tumour samples
################################################################################

#suffixs = c("149", "70s", "50s")
#suffixs = c("149", "70s")
suffixs = c("149")
for(suffix1 in suffixs){
  if(suffix1 == "70s"){
      patient_info_tmp <-
        patient_info[rho_reseg >= 0.7 & rho_reseg <= 0.75]
      }
  if(suffix1 == "50s"){
      patient_info_tmp <-
        patient_info[rho_reseg >= 0.52 & rho_reseg <= 0.57]
      }
  if(suffix1 == "149"){
      patient_info_tmp <-
        patient_info[Patient=="LTX149"]
      }

  patient_ids <- unique(patient_info_tmp$Patient)
  dt_list <- list()
  for(patient_id in patient_ids) {
    path_patient <- paste0(tx_dir, patient_id, "/",sep="")
    
    # get matched normal info
    f_nm <- paste0(path_patient,"/Allelecounts/N/",
                   patient_id,".N.SNPs.CpGs.all.sorted.RData")
    load(f_nm)
    gr_normal$mq <- NULL
    dt <- data.table(data.frame(gr_normal, stringAsFactors = FALSE))
    rm(gr_normal, f_nm)
    
    #Get SNP in ref:alt format 
    dt[, SNP := ifelse(!is.na(ref), paste(ref, alt, sep = ""), NA)]
    dt[, offset := ifelse(width == 4, 1, 0)]
    
    # Annotate SNP type
    dt[, "type" := factor(ifelse(BAF < 0.15,"Homozygous",
                        ifelse(BAF > 0.85,"Homozygous",
                        ifelse(BAF >= 0.15 & BAF <= 0.85, "Heterozygous",
                               NA))), 
                        levels = c("Homozygous", "Heterozygous"))]

    # Extract columns of interest
    colnames(dt)[which(colnames(dt)=="seqnames")] <- "chromosome"
    dt_norm <- dt[width > 1 & !is.na(total_counts_m) & total_counts_m >= 10, ]
    cols <- c("width", "chromosome", "start", "end", 
              "SNP", "offset", "POS", "BAF", "type")
    dt_norm <- dt_norm[, .SD, .SDcols = cols]
    rm(dt)
    
    # get tumour sample ids
    sample_ids <- patient_info_tmp$Sample[patient_info_tmp$Patient == patient_id]
    
    # get tumour sample info
    for (sample_id in sample_ids) {
      print(patient_id)
      print(sample_id)
      
      # load the tumour data
      f_nm <-
        paste0(path_patient, "/Allelecounts/", sample_id, "/",
              patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
      load(f_nm)
      gr_tumour$mq <- NULL
      dt <- data.table(data.frame(gr_tumour, stringAsFactors=FALSE))

      # rename BAF column
      dt[, "tumour_BAF" := BAF]

      # subset to CpGs
      colnames(dt)[which(colnames(dt) == "seqnames")] <- "chromosome"
      cols <- c("width", "chromosome", "start", "end", "POS", "tumour_BAF",
                "Tr", "Af", "Tf", "Ar",
                "CGf", "CGr", "TGf", "CAr")
      dt <- dt[, .SD, .SDcols = cols]
      dt <- dt[width > 1, ]
      rm(gr_tumour)
      
      # assign normal BAF
      dt <-
        merge(dt, dt_norm, by=c("chromosome", "start", "end", "width", "POS"), all = FALSE)

      # Counts methylated for allele-specific methylation i.e. excluding CG destroying SNPs.
      dt[, M_CAMDAC := 
                ifelse(is.na(SNP), CGf + CGr,
                ifelse(start + offset == POS & 
                          grepl(pattern = "([CT][TC])",SNP), CGr,
                ifelse(end - offset == POS & 
                          grepl(pattern = "([GA][AG])",SNP), CGf,
                        CGf + CGr)))]
      dt[, M_NAIVE := CGf + CGr] #i.e. Cf and Gr

      # Get the total methylation counts --> CG-destroying SNP allele does not contribute. 
      # Will be NA if the CG-destroying allele at CG/CCGG is the only one present
      dt[, "UM_CAMDAC" := 
                ifelse(is.na(SNP),
                  TGf + CAr,
                ifelse(start + offset == POS &
                  grepl(pattern = "([CT][TC])", SNP),
                    CAr,
                ifelse(end - offset == POS &
                  grepl(pattern = "([GA][AG])", SNP),
                    TGf,
                    TGf + CAr)))]
      dt[, "UM_NAIVE" :=
                ifelse(is.na(SNP),
                  TGf + CAr, #i.e. Cf and Ar
                ifelse(
                  (start + offset == POS &
                    grepl(pattern = "(CA|CG|AC|GC)", SNP)),
                  TGf + CAr + Ar, # Ar = AAr 
                ifelse(
                  (end - offset == POS &
                    grepl(pattern = "(TG|CG|GT|GC)", SNP)),
                  TGf + CAr + Tf, # Tf = TTf 
                ifelse(start + offset == POS &
                  grepl(pattern = "(CT|TC)", SNP),
                  TGf + CAr + Tr, # Tr representing TAr
                ifelse(end - offset == POS &
                  grepl(pattern = "(GA|AG)", SNP),
                  TGf + CAr + Af, # Af representing TAf
                NA)))))]

      # Get the counts unmethylated
      dt[, "total_counts_m_CAMDAC" := UM_CAMDAC + M_CAMDAC]
      dt[, "total_counts_m_NAIVE" := UM_NAIVE + M_NAIVE]

      # Remove low coverage loci
      dt <- dt[total_counts_m_CAMDAC >= 3 & total_counts_m_NAIVE >= 3,]

      # Get the counts unmethylated
      dt[, m_CAMDAC := M_CAMDAC / total_counts_m_CAMDAC]
      dt[, m_NAIVE := M_NAIVE / total_counts_m_NAIVE]
      dt[, `:=`(offset=NULL, M_CAMDAC = NULL, UM_CAMDAC = NULL, M_NAIVE = NULL, UM_NAIVE = NULL)]

      # Set color flag for SNP type
      dt[, "SNP_type" := ifelse(is.na(SNP) | type == "Homozygous", "None", 
        ifelse(grepl("(CT)|(TC)|(GA)|(AG)", SNP), "CT", "CAorCG"))]
      cols <- c("CT" = "darkblue", "CAorCG" = "gold3", "None" = "grey50")
    
      # Load the copy number data
      load(paste0("~/rdeocd/txMeth/",patient_id,"/Copy_number/",sample_id,
          "/re_segmented/",patient_id,".",sample_id,".ascat.output.resegmented.RData"))
      seg <- data.table(ascat.output$segments[, c("chr", "startpos", "endpos", "nMajor", "nMinor")])
      seg[, "chromosome" := as.character(chr)]
      seg[, "start" := as.integer(startpos)] 
      seg[, "end" := as.integer(endpos)]
      rm(ascat.output)
      
      # Annotate copy number onto tumour data.table
      dt[, chromosome := as.character(substr(chromosome, 4, 5))]
      setkey(dt, chromosome, start, end)
      setkey(seg, chromosome, start, end)
      dt_seg <- foverlaps(dt, seg, nomatch = 0)
      dt_seg[, `:=`(start = i.start, end = i.end)]
      dt_seg[, `:=`(i.start = NULL, i.end = NULL, chr = NULL, startpos = NULL, endpos = NULL)]
            
      # Store sample data
      dt_list[[length(dt_list)+1]] <- dt_seg
      rm(dt_seg, dt, seg)
    }
  }
  # Store group data
  assign(paste0("dt_list_", suffix1), dt_list)
}

# Save outputs
save(dt_list_149, file = file.path("outputs", "combined_tumours.Rdata"))
rm(dt_norm, dt_list)

# Plot the comparison between NAIVE and CAMDAC methylation rates
for (suffix1 in suffixs) {
  dt_list <- mget(paste0("dt_list_", suffix1))[[1]]
  dt_seg <- rbindlist(dt_list)
  dt_seg[, nTotal := nMajor + nMinor]
    
  suffix <- suffix1
  tmp2 <- dt_seg

  # create empty plotting square
  empty <- ggplot() + geom_point(aes(1, 1), colour = "white")+
         theme(axis.ticks = element_blank(), panel.background = element_blank(), 
               axis.text.x = element_blank(), axis.text.y = element_blank(),           
               axis.title.x = element_blank(), axis.title.y = element_blank())

  for(i in 1:3){
    # Subset by tumour copy number
    tmp <- tmp2[nTotal == i]
    nT <- i

    # Plot all CpGs
    d_top <- ggplot(tmp) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
    scatter <- ggplot(tmp) + theme_classic() +
              geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
                alpha = 0.1, shape = 16, size = 1) +
              scale_color_manual(values = cols) +
              theme(legend.position = "none")
    d_right <- ggplot(tmp) + theme_classic() + 
              geom_histogram(
                aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
                fill=NA,  binwidth=0.04, alpha=0.1) +
              scale_color_manual(values = cols) +
              coord_flip() +
              theme(axis.text.x = element_text(angle = 270),
                    legend.position = "none")
    gd <-
      grid.arrange(d_top, empty, scatter, d_right, 
                  ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_", nT,"_",suffix,".png")),
             plot = gd, device = "png", width=5, height=5)
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_", nT, "_", suffix, ".pdf")),
           plot = gd, device = "pdf", width=5, height=5)
    
    # subset to all heterozygous SNPs
    tmp <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7]
    d_top <- ggplot(tmp) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
    scatter <- ggplot(tmp) + theme_classic() +
              geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
                alpha = 0.1, shape = 16, size = 1) +
              scale_color_manual(values = cols) +
              theme(legend.position = "none")
    d_right <- ggplot(tmp) + theme_classic() + 
              geom_histogram(
                aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
                fill=NA,  binwidth=0.04, alpha=0.1) +
              scale_color_manual(values = cols) +
              coord_flip() +
              theme(axis.text.x = element_text(angle = 270),
                    legend.position = "none")
    gd <-
      grid.arrange(d_top, empty, scatter, d_right, 
                  ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_", nT, "_", suffix, ".png")),
           plot = gd, device = "png", width=5, height=5)
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_", nT, "_", suffix, ".pdf")),
           plot = gd, device = "pdf", width=5, height=5)
    
    # subset to C>T nd G>A SNPs
    tmp3 <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7 &
      grepl("(CT)|(TC)|(GA)|(AG)", SNP)]
    d_top <- ggplot(tmp3) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
    scatter <- ggplot(tmp3) + theme_classic() +
              geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
                alpha = 0.1, shape = 16, size = 1) +
              scale_color_manual(values = cols) +
              theme(legend.position = "none")
    d_right <- ggplot(tmp3) + theme_classic() + 
              geom_histogram(
                aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
                fill=NA,  binwidth=0.04, alpha=0.1) +
              scale_color_manual(values = cols) +
              coord_flip() +
              theme(axis.text.x = element_text(angle = 270),
                    legend.position = "none")
    gd <-
      grid.arrange(d_top, empty, scatter, d_right, 
                  ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_CT_", nT, "_", suffix, ".png")),              
           plot = gd, device = "png", width=5, height=5)
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_CT_", nT, "_", suffix, ".pdf")),                        
           plot = gd, device = "pdf", width=5, height=5)

    # subset to C>A nd C>G SNPs
    tmp3 <- tmp[!is.na(SNP) & BAF >= 0.3 & BAF <= 0.7 &
      grepl("(CA)|(AC)|(GC)|(CG)", SNP)]
    d_top <- ggplot(tmp3) + theme_classic() +
         geom_histogram(
          aes(x = m_NAIVE, y = ..count.., col = SNP_type),
          fill = NA, binwidth = 0.04, alpha = 0.1) +
         scale_color_manual(values = cols) +
         theme(axis.text.x = element_text(angle = 270),
              legend.position = "top")
    scatter <- ggplot(tmp3) + theme_classic() +
              geom_point(aes(m_NAIVE, m_CAMDAC, col = SNP_type),
                alpha = 0.1, shape = 16, size = 1) +
              scale_color_manual(values = cols) +
              theme(legend.position = "none")
    d_right <- ggplot(tmp3) + theme_classic() + 
              geom_histogram(
                aes(x = m_CAMDAC, y = ..count.., col = SNP_type),
                fill=NA,  binwidth=0.04, alpha=0.1) +
              scale_color_manual(values = cols) +
              coord_flip() +
              theme(axis.text.x = element_text(angle = 270),
                    legend.position = "none")
    gd <-
      grid.arrange(d_top, empty, scatter, d_right, 
                  ncol = 2, nrow = 2, widths = c(2, 1), heights = c(1.5, 2))
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_CAorCG_", nT, "_", suffix, ".png")),              
           plot = gd, device = "png", width=5, height=5)
    ggsave(filename = file.path("outputs",
      paste0("tumour_methylation_SNP_CAorCG_", nT, "_", suffix, ".pdf")),                        
           plot = gd, device = "pdf", width=5, height=5)

    # clean-up
    rm(d_top, d_right, gd, scatter, tmp3)
  }
}

# clean-up
rm(dt_seg, tmp, tmp2)
rm(list = ls(pattern = "dt_list"))

# END