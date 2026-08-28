### Most updated RRBS imprinting analysis code adapted from Lili's code and using Lili's input files 
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
df_imprinted_genes <- as.data.table(read_tsv("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/Imprinted_genes_2021.txt"))

# get chromosomal coordinates from GRCH37 to hg19
df_imprinted_genes[, chromosome := unlist(tstrsplit(Genome_Location, ":", keep=1))]
df_imprinted_genes[, start := tstrsplit(unlist(tstrsplit(Genome_Location, ":", keep=2)), "-", keep=1)]
df_imprinted_genes[, end := tstrsplit(unlist(tstrsplit(Genome_Location, ":", keep=2)), "-", keep=2)]
df_imprinted_genes[, start := as.numeric(start)]
df_imprinted_genes[, end := as.numeric(end)]
# df_imprinted_genes[, start := as.numeric(start)-5000]
# df_imprinted_genes[, end := as.numeric(end)+5000]
df_imprinted_genes[, Genome_Location_tmp := NULL]

# correct H19 coords
df_imprinted_genes[Gene=="H19", start := 2018500]
df_imprinted_genes[Gene=="H19", end := 2022500]

# set genome coordinates key
setkeyv(df_imprinted_genes, cols=c("chromosome", "start", "end"))

# remove the imprinted loci covering cytogenic bands for now
df_imprinted_genes <- df_imprinted_genes[!is.na(end)]

# only keep clear Imprinted genes either maternal, paternal or non-imprinted
df_imprinted_genes <- df_imprinted_genes[!Expressed_Allele%in%c("Biallelic", "Unknown")]
df_imprinted_genes <- df_imprinted_genes[Status=="Imprinted"]
genes <- df_imprinted_genes$Gene
write_csv(df_imprinted_genes, "Imprinted_lili_hg19.csv")

# convert hg19 to hg 38 with liftover
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
CN_good <- readLines("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/Good_samples_59.txt")
metadata <- read.csv("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/CAMDAC_metadata_MSS.csv")
meta_rrbs <- read.csv("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/rrbs_cn_summary.csv")
metadata <- inner_join(metadata, meta_rrbs, by=c("PT_ID"="sample_id"))
metadata <- metadata %>% filter(PT_ID %in% CN_good)
asm_paths_rrbs <- Sys.glob(paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/final_version/ASM_RRBS/", metadata$PT_ID, "/N/allele_specific_methylation_rates_N.RData"))
bulk_normal_path <- Sys.glob(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/Batch*/bam_final/", metadata$PT_ID, "/Methylation/N1/dt_normal_m.RData"))

# get consectuvive DMP helper function
fct <- function(x,sgn,min) {
  flag <- c(FALSE, diff(sign(sgn)) > 0)
  y <- rle(x==TRUE & flag==FALSE)
  z<-as.numeric(suppressWarnings(
    max(y$lengths[y$lengths >= min & y$values == TRUE]))) 
  if(!is.finite(z)){z <- as.numeric(NA)}
  return(z)
} 

# get imprinted bins
dt_IM_regions_list <- list()
for(i in 1:nrow(metadata)){
  load(asm_paths_rrbs[i]) 
    
  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, ov)
  df_phased_methylation <- df_phased_methylation[!is.na(CHR)&!is.na(Gene)]
    
  # annotate imprinted gene loci
  #df_phased_methylation[, as.character(c("i.start", "i.end")) := NULL]
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
table(dt_IM_regions$IM_m_n_phased) # 159 bins with imprinting; there are overlaps
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
setkeyv(dt_IM_regions_anno, cols=c("chrom", "start", "end"))

# get normal methylation at imprinted bins
dt_IM_regions_normal_list <- list()

for(i in 1:nrow(metadata)){
  load(asm_paths_rrbs[i])

  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, dt_IM_regions_anno)
  df_phased_methylation <- df_phased_methylation[!is.na(CHR)&!is.na(Gene)]
  
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
  df_phased_methylation[, patient := metadata$PT_ID[i]]
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
  load(bulk_normal_path[i])
  df_bulk <- data.table(data.frame(
    dt_normal_m))[, .("chrom"=substr(as.character(CHR),4,5), 
                      start, end, m_n, m_n_low, m_n_high)]
  rm(dt_normal_m)
  
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
      geom_text(data=tmp[1,],aes(x=unique(as.numeric(POS_genotype)),y=-0.1,label = "*"),col="firebrick",size=10)+
      theme(axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black")) + 
      scale_color_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))+
      scale_fill_manual(name="",values=c("ref"="gold2", "alt"="mediumpurple1"))
    
    p1 <- p + geom_point(data=df_anno, aes(x=start,y=m_n),col="grey10")+
      geom_ribbon(data=df_anno, aes(x=start,ymin=m_n_low,ymax=m_n_high), fill="grey10", alpha=.2)+
      geom_point(data=tmp, aes(x=CpG_start,y=m_n_1,col="ref"))+
      geom_point(data=tmp, aes(x=CpG_start,y=m_n_2,col="alt"))
    
    ggsave(filename =  paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/final_version/ASM_RRBS/", metadata$PT_ID[i], "/N/", gene, "_", ID,"_plot.pdf"), 
           plot=p1, device = "pdf", units = "in", width = 5, height = 3)
    rm(p1, p, df_anno, tmp, coords, p1)
  }
  rm(df_bulk, df_phased_methylation)
} 
dt_IM_normals <- rbindlist(dt_IM_regions_normal_list)
save(dt_IM_normals, dt_IM_regions_anno, file="ASM/normal_imprinting_20240627.RData")

library(purrr)
merged_df <- dt_IM_regions_normal_list %>%
  map(~ select(.x, Gene, patient)) %>%
  bind_rows()

# create a heatmap for imprinting 
df_wide <- merged_df %>% distinct() %>%
  mutate(value = 1) %>%  # Create a column to fill with 1s
  pivot_wider(names_from = patient, values_from = value, values_fill = list(value = 0)) %>%
  mutate(across(-Gene, ~ ifelse(is.na(.), 0, .)))

mat <- as.matrix(df_wide[,-1])
rownames(mat) <- df_wide$Gene

row_sums <- rowSums(mat)
# order_index <- order(-row_sums)
order_index <- order(rownames(mat))
mat <- mat[order_index, ]

library(pheatmap)
p <- pheatmap(mat, 
         color = c("white", "red"), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend_breaks = c(0,1),
         legend_labels = c("unknown", "imprinting"))

save_pheatmap_pdf <- function(x, filename, width=12, height=4) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

save_pheatmap_pdf(p, "ASM/rrbs_normal_consec_lili_updated.pdf")



# imprinting in tumors -----------
load("ASM/normal_imprinting_20240627.RData")
dt_IM_tumours_list <- list()
asm_paths_rrbs <- Sys.glob(paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/final_version/ASM_RRBS/", metadata$PT_ID, "/T/allele_specific_methylation_rates_T.RData"))

for(1 in 1:nrow(metadata)){
  print(i)
  load(asm_paths_rrbs[i]) 
  
  # annotate imprinted gene loci
  setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
  df_phased_methylation <- foverlaps(df_phased_methylation, dt_IM_regions_anno)
  df_phased_methylation <- df_phased_methylation[(!is.na(m_t_1)|!is.na(m_t_2))&!is.na(Gene)]
  
  # filter out CpGs with insufficient counts and/or non-CG SNPs
  df_phased_methylation <- df_phased_methylation[is.na(label)]
  
  # remove duplicate colnames
  colnames(df_phased_methylation) <- gsub("i\\.", "CpG_", colnames(df_phased_methylation))
  df_phased_methylation <- df_phased_methylation[!is.na(Gene)][order(id, Gene, CpG_start, CpG_end)]
  
  # filter out bins IM and validated by SNP in normal
  tmp <- dt_IM_normals[patient==metadata$PT_ID[i]&IM_m_n_phased=="IM"]
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
  tmp <- df_phased_methylation[, .SD, .SDcols=c("m_diff_t12", "IM_t_average", "id")]
  tmp <- tmp[,.(consec_IM = fct(x = IM_t_average=="IM", sgn=m_diff_t12, min = 2)), by=.(id)]
  df_phased_methylation$consec_IM_m_t_average <-
      tmp[match(df_phased_methylation$id, id), .SD, .SDcols=c("consec_IM")]
  rm(tmp)

  # Annotate imprinted genes according to consecutive bulk IM CpGs
  df_phased_methylation[, IM_m_t_average := ifelse(is.na(m_t_1)&is.na(m_t_2), "no_info", 
                                            ifelse(!is.na(consec_IM_m_t_average) &
                                                    consec_IM_m_t_average >= 2, "IM","not_IM"))]
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
  df_phased_methylation[, patient := metadata$PT_ID[i]]
  df_phased_methylation[, LOH := ifelse((n1==0|n2==0)&(n2_n!=0&n1_n!=0), 1, 0)]
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
  # dt_info[, c("mean_m_t","mean_m_t_phased") := 
  #           list(ifelse(IM_m_t_phased=="no_info", as.character(NA),
  #                        ifelse(m_t>0.5, "M", "UM")),
  #               ifelse(IM_m_t_phased=="no_info", as.character(NA),
  #               ifelse((n2<0.75 & m_t_1>0.5) | (n1<0.75 & m_t_2>0.5) | 
  #                      (n2>=0.75&n1>=0.75&(n1*m_t_1+n2*m_t_2)/(n1+n2)>0.5), "M",
  #               ifelse((n2<0.75 & m_t_1<=0.5) | (n1<0.75 & m_t_2<=0.5)| 
  #                      (n2>=0.75&n1>=0.75&(n1*m_t_1+n2*m_t_2)/(n1+n2)<=0.5), "UM",
  #                       as.character(NA)))))]
    dt_info[CpGs_t < 3, IM_m_t_phased := "no_info"]
    dt_info[CpGs_t < 3, IM_m_t_average := "no_info"]
    dt_info[, flag := ifelse(IM_m_t_phased == "no_info", "no_info", 
                      ifelse(LOH == TRUE, "LOH", 
                      ifelse(IM_m_t_phased == "IM", "Imprinted_phased",
                      ifelse(m_diff_t12 >= 0.5, "Imprinted_maybe", 
                      ifelse(m_diff_t12 <= 0.15 & IM_m_t_phased == "not_IM", "LOI", 
                      ifelse(m_diff_t12 < 0.5, "partial_LOI", "investigate"))))))]

    # Load binned info
    file_path=paste0(path_patient, "/Methylation/", sample_id, "/re_segmented/")
    file = paste0(file_path, "purified_tumour_DMR_calls.RData")
    if(!file.exists(file)){
      file_path=paste0(path_patient, "/Methylation/", sample_id, "/raw/")
      file = paste0(file_path, "purified_tumour_DMR_calls.RData")
    }
    load(file)
    df_regions <- data.table(df_annotated_regions) ; rm(df_annotated_regions)
    dt_info[, c("m_n_bin", "m_t_bin") := df_regions[match(dt_info$id, id),  .(m_n, m_t)]]
    
    dt_IM_tumours_list[[length(dt_IM_tumours_list)+1]] <- dt_info
    rm(dt_info, df_regions, file, file_path)
} ; rm(df_phased_methylation)

dt_IM_tumours <- rbindlist(dt_IM_tumours_list)
save(dt_IM_tumours, file = "~/rdeocd/txMeth/imprinting_summary_tumours_Oct2021_v3.RData")
rm(dt_IM_tumours_list)


# get number of phased CpGs to SNPs
PT_ID <- c("LORC-2", "LORC-60", "LORC-62", "LORC-68", "YORC-25", "YORC-41")
asm_paths_n <- Sys.glob(paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/", PT_ID, "/N/allele_specific_methylation_rates_N.RData"))
asm_paths_t_v2 <- Sys.glob(paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/", PT_ID, "/T/allele_specific_methylation_rates_T.RData"))
asm_paths_t <- Sys.glob(paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/", PT_ID, "/T/allele_specific_methylation_rates_T_v1.RData"))

# normal
for (f in asm_paths_n){
  load(f)
  print(nrow(df_phased_methylation))
  rm(df_phased_methylation)
}
# tumor v2 as named in lili's original asm code
for (f in asm_paths_t_v2){
  load(f)
  print(nrow(df_phased_methylation))
  rm(df_phased_methylation)
}

# tumor newest version as named in lili's original asm code last section
for (f in asm_paths_t){
  load(f)
  print(nrow(df_phased_methylation))
  rm(df_phased_methylation)
}
