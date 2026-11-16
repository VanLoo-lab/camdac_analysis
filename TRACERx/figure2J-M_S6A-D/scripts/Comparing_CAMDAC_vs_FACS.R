#' @title Compare FACS-purified and CAMDAC methylation rates
#' @author Lili Larose Cadieux
#' @description In this code, we compare:
#' * CAMDAC bulk tumour- adjacent normal vs. FACS purified
#'    tumour-infiltrating normals
#' * CAMDAC purified tumour vs. FACS purified tumour
#' * CAMDAC bulk tumour vs. FACS purified tumour

# Load libraries
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
x <- c("stringr", "GenomicRanges", "scales", "ggplot2", "gridExtra",
       "data.table", "readr", "fst")
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

# set parent data directories
parent_dir <- "~/rdeocd/FACS_sorted_txMeth/data/"
parent_dir2 <- "~/rdeocd/txMeth/"
work_dir <- "~/data/CAMDAC_paper/figure2J-M_S6A-D"

# load patient metadata
f_nm <- "NuGEN_summary_report_FACS_sorted_samples.txt"
patient_info <-
  fread(file = paste0(parent_dir, f_nm), header = TRUE, sep = "\t")

# remove "," separator from large numbers (excel formatting artifact)
# remove comma separators (excel formatting artifact) and convert
# numeric columns to proper type
numeric_cols <- colnames(patient_info)[7:ncol(patient_info)]
patient_info <- patient_info[,
  # Remove commas from all columns
  lapply(.SD, function(x) gsub(pattern = ",", replacement = "", x)),
][,
  # Convert specified columns to numeric
  (numeric_cols) := lapply(.SD, as.numeric),
  .SDcols = numeric_cols
]

# annotate FACS population (1 = Normal, 2 = tumor)
patient_info[(nchar(as.character(Patient)) == 6L & substr(as.character(Patient), 4, 4) == "0"), pop_id := substr(Sample_Name, 7, 10)]
patient_info[!(nchar(as.character(Patient)) == 6L & substr(as.character(Patient), 4, 4) == "0"), pop_id := substr(Sample_Name, 8, 11)]

# load clinical patient_info
f_nm <- file.path(parent_dir2,
               "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt")
patient_clin <- fread(f_nm, header = TRUE, sep = "\t")
rm(f_nm)

# remove lymph nodes
patient_clin <- patient_clin[!grepl("LN", Sample), ]

# annotate CRUKID patient ID
patient_info[, CRUKID := patient_clin[match(patient_info$Patient, Patient), CRUKID]]
patient_info[patient_clin, "CRUKID2" := i.CRUKID, on = .(Patient)]

# annotate sample ID
patient_info[, LTXID := paste0(patient_info$Patient, "_", patient_info$Region)]

# annotate purity and ploidy
patient_info[, rho := patient_clin[match(patient_info$LTXID, LTX_ID), rho_reseg]]
patient_info[, psi := patient_clin[match(patient_info$LTXID, LTX_ID), psi_reseg]]

################################################################################
# Plot FACS purified vs. CAMDAC bulk normals methylation rate
################################################################################

# create object to store interesting results
dt_list <- list()
dt_list2 <- list()
dt_list3 <- list()

# loop over each normal sample
for (sample_id in patient_info[, sort(Sample_Name)]){
  patient_id <- patient_info[Sample_Name == sample_id, Patient]
  region_id <- patient_info[Sample_Name == sample_id, Region]
  pop_id <- paste0("pop", patient_info[Sample_Name == sample_id, Pop])

  if (pop_id == "pop1") {
    print(sample_id)

    # create dirs
    path <- paste0(parent_dir, patient_id, "/")

    # get file name prefix
    path_output <- paste0(path, "Methylation/", sample_id, "/")

    # set filenames and filepaths
    output_file <- paste0(path_output, "purified_tumour_DMP_calls.RData")

    # load file processed methylation info
    load(output_file)
    dt <- data.table(df_mt)[, .SD, .SDcols = c("chromosome", "start", "end",
                                               "m_n", "m_n_low", "m_n_high",
                                               "m_b", "m_b_low", "m_b_high",
                                               "DMP_b", "DMP_b_corr")]
    suppressWarnings(rm(gr_sample_v_normal_m, df_mt))
    colnames(dt) <- gsub("m_b", "m_f", colnames(dt))
    colnames(dt) <- gsub("m_n", "m_b", colnames(dt))

    # class CpGs
    dt[, flag := ifelse(m_f_high < 0.2, "UM",
                 ifelse(m_f_low > 0.8, "FM",
                 ifelse(m_f_low > 0.2 & m_f_high < 0.8, "IM",
                          "NA")))]
    # class CpGs
    dt[, flag2 := ifelse(m_b_high < 0.2, "UM",
                  ifelse(m_b_low > 0.8, "FM",
                  ifelse(m_b_low > 0.2 & m_b_high < 0.8, "IM",
                           "NA")))]
    # class CpGs
    dt[, flag3 := ifelse(m_b_high < 0.2 | m_f_high < 0.2, "UM",
                  ifelse(m_b_low > 0.8 | m_f_low > 0.8, "FM",
                  ifelse((m_b_low > 0.2 & m_b_high < 0.8) |
                         (m_f_low > 0.2 & m_f_high < 0.8), "IM",
                          "NA")))]

    # load cell type markers
    dt_cell_types <-
      read.table(file = "~/rdeocd/pipeline_files/cell_type_markers.txt",
                 sep = "\t", header = TRUE)
    dt_cell_types <- data.table(dt_cell_types)[RRBS > 0,]
    dt_cell_types[, m_b_n := as.numeric(NA)]
    dt_cell_types[, m_f_n := as.numeric(NA)]
      
    # Load Bulk vs.FACS DMRs
    f_name <- paste0("~/rdeocd/FACS_sorted_txMeth/data/", patient_id,
                       "/Methylation/", sample_id, "/bulk_tumour_DMR_calls.RData")
    load(f_name)
      
    # select relevant columns
    df_regions <-
      df_annotated_regions[, .SD,
      .SDcols = c("chromosome", "start", "end", "m_n_DMPs", "m_b_DMPs", "prob_b_DMPs",
                  "CGI", "gene_ids","gene_prom_ids","enhancer_ids",
                  "DMP_b", "consec_DMPs", "DMR")]
    rm(f_name, df_annotated_regions)

    # convert to numeric
    df_regions[, as.character(c("m_n_DMPs", "m_b_DMPs", "prob_b_DMPs")) :=
                   lapply(.SD, as.numeric),
                 .SDcols = c("m_n_DMPs", "m_b_DMPs", "prob_b_DMPs")]

    # rename appropriately
    colnames(df_regions) <- gsub("m_b", "m_f", colnames(df_regions))
    colnames(df_regions) <- gsub("m_n", "m_b", colnames(df_regions))
    df_regions[, m_diff_DMPs := m_b_DMPs - m_f_DMPs]
    df_regions[, alphaa := ifelse(gene_prom_ids!="", "yes", "no")]
    # df_regions[, alphaa := ifelse(enhancer_ids!="", "yes", "no")]
    # df_regions[, alphaa := ifelse(CGI!="NA", "yes", "no")]

    ggplot(df_regions[!is.na(DMR),]) + theme_minimal() +
      geom_point(data=df_regions[!is.na(DMR)&DMP_b=="hyper"&consec_DMPs>5,],
                 aes(x=m_diff_DMPs, y=-log10(prob_b_DMPs), alpha=alphaa), col="red")+
      geom_point(data=df_regions[!is.na(DMR)&DMP_b=="hypo"&consec_DMPs>5,],
                 aes(x=m_diff_DMPs, y=-log10(prob_b_DMPs), alpha=alphaa), col="blue")+
      geom_point(data=df_regions[!is.na(DMR)&DMP_b=="mixed"&consec_DMPs>9,],
                 aes(x=m_diff_DMPs, y=-log10(prob_b_DMPs), alpha=alphaa), col="grey65")+
      scale_alpha_manual(values=c("yes"=1, "no"=0.1))

    for (gene in unique(dt_cell_types$Marker)) {
      y <- paste0("(;", gene, ",)|(^", gene, ",)", collapse = "|")
      tmp <- df_regions[grepl(y, gene_prom_ids) & !is.na(consec_DMPs), ]
      if (nrow(tmp) > 0) {
        DMR_flag <- tmp[, sum(!is.na(DMR), na.rm = TRUE) > 0]
        dt_cell_types[Marker == gene, DMR := DMR_flag]
        if (DMR_flag == TRUE) {
          hyper_flag <- tmp[, ifelse(sum(DMP_b == "hyper") > 0 &
                                       sum(DMP_b == "hypo") > 0, "mixed",
                                ifelse(sum(DMP_b == "hyper") > 0, "hyper", "hypo"))]
          if (nrow(tmp) == 1){
            dt_cell_types[Marker == gene, m_b_n := tmp[, m_b_DMPs]]
            dt_cell_types[Marker == gene, m_f_n := tmp[, m_f_DMPs]]
          }
          if (nrow(tmp) > 1 & hyper_flag == "mixed") {
            idx <- which.max(tmp$consec_DMPs)
            dt_cell_types[Marker == gene, m_b_n := tmp[idx, m_b_DMPs]]
            dt_cell_types[Marker == gene, m_f_n := tmp[idx, m_f_DMPs]]
          }
          if (nrow(tmp) > 1 & hyper_flag != "mixed") {
            sumss <- abs(sum(tmp$consec_DMPs, na.rm = TRUE))
            dt_cell_types[Marker == gene, m_b_n :=
                            tmp[, sum(m_b_DMPs * consec_DMPs, na.rm = TRUE) /
                            sumss]]
            dt_cell_types[Marker == gene, m_f_n :=
                            tmp[, sum(m_f_DMPs * consec_DMPs, na.rm = TRUE) /
                            sumss]]
          }
        }
      }
    }
    dt_cell_types[, m_diff := m_b_n - m_f_n]
    dt_cell_types[, patient := patient_id]
    dt_cell_types[, sample := sample_id]
    dt_list3[[length(dt_list3) + 1]] <- dt_cell_types[m_diff != 0]
    # # overlap DMRs with IM CpG data
    # tmp <- dt[(flag=="IM" | flag2 =="IM"), .SD,
    #           .SDcols=c("chromosome", "start", "end","flag","flag2")]
    # setkeyv(df_regions, c("chromosome", "start", "end"))
    # tmp2 <- foverlaps(tmp, df_regions)
    # tmp2 <- tmp2[!is.na(m_diff_DMPs),]
    #
    # # select all hyper
    # tmp2_hyper <- tmp2[(DMP_b=="hyper"),]
    # tmp2_hyper <- tmp2_hyper[, .(m_diff_DMPs=m_diff_DMPs[1]),
    #                          by=.(chromosome, start, end,
    #                               gene_ids, gene_prom_ids, CGI,
    #                               DMR)]
    # v <- tmp2_hyper[, tstrsplit(gene_prom_ids, ";", fixed=TRUE)]
    # w <- tmp2_hyper[, tstrsplit(gene_prom_ids, ";", fixed=TRUE)]
    # w <- c(v, w)
    # x <- data.table(w=unique(na.omit(unlist(w))))
    # z <- x[, y := tstrsplit(w, ",", fixed=TRUE, keep=2)]$y
    # cat(z, sep=" ")
    # z <- x[, y := tstrsplit(w, ",", fixed=TRUE, keep=1)]$y
    # cat(z, sep=" ")
    # rm(tmp, tmp2, v, w, x, y, z)

    # set the image output file prefix
    f_prefix <- paste0(work_dir, "/outputs/", sample_id, "_FACS_normal_vs_CAMDAC_")

    # plot the methylation histogram
    df_cor <- data.frame(
      x = 0.85,
      y = 1E6,
      tex = paste0("cor=", dt[, round2(cor(m_b, m_f), digits = 3)])
      )
    p <- ggplot(dt, aes(x=m_b, y = after_stat(count))) + theme_classic() +
      geom_histogram(aes(col = "Bulk", fill = "Bulk"), alpha = 0.1, bins = 25) +
      geom_histogram(aes(col = "FACS", fill = "FACS", x = m_f),
        alpha = 0.1, bins = 25) + 
      ylab("counts") + ggtitle(sample_id) +
      geom_text(data = df_cor, aes(x = x, y = y, label = tex)) +
      scale_color_manual(name = "platform",
                         values = c("Bulk" = "forestgreen", "FACS" = "purple3")) +
      scale_fill_manual(name = "platform",
                        values = c("Bulk" = "forestgreen", "FACS" = "purple3"))
    ggsave(plot = p, filename = paste0(f_prefix, "_bulk_normal_hist.pdf"),
           device = pdf, width = 4, height = 3, unit = "in")

    # set boundaries for zoomed-in figures
    q1 <- dt[flag3 != "NA", quantile(m_b - m_f, 0.025)]
    q2 <- dt[flag3 != "NA", quantile(m_b - m_f, 0.975)]

    # get correlations
    cor_val <- paste0("cor=", dt[flag3 != "NA", round2(cor(m_b, m_f), digits = 3)])
    
    p2 <- ggplot(data = dt[flag3 != "NA"]) +  theme_classic() +
      geom_histogram(aes(x = m_b - m_f, y = after_stat(count),
                         fill = factor(flag3, levels = c("FM", "IM",  "UM"), ordered = TRUE),
                         alpha = factor(flag3, levels = c("FM", "IM", "UM"), ordered = TRUE)),
                     bins = 100, position = "identity", col = NA) +
      ggtitle(sample_id) +
      #geom_vline(xintercept = c(q1, q2), linetype="dashed", col="black") +
      geom_hline(yintercept = 20000, linetype="dashed", col="black") +
      geom_text(data=data.frame(x = 0.75, y = 2E2, tex = cor_val),
                  aes(x=x, y=y, label=tex)) +
      theme(panel.grid.major = element_line(size=.1, color="black")) +
      scale_x_continuous(name = "mn - mf,n", breaks=seq(-1, 1, 0.25), limits = c(-1, 1)) +
      scale_y_continuous(name = "counts", breaks=seq(0, 8E5, 2E5), 
                           limits = c(0, 800000), labels = paste0(seq(0, 800, 200), "k")) +
        scale_fill_manual(name="state", values=c("UM" = "blue", "IM" = "darkturquoise", "FM" = "gold2")) +
      scale_alpha_manual(name="state", values=c("UM" = 0.35, "IM" = 0.75, "FM" = 0.35))
    ggsave(plot = p2, filename = paste0(f_prefix, "_vs_normal_hist_diff_all.pdf"), device = "pdf",
             width=4, height = 3, unit="in")
      
    p2_zoom <-  ggplot(data = dt[flag3 != "NA"]) +  theme_classic() + 
        geom_histogram(aes(x = m_b - m_f, y = after_stat(count), 
                           fill = factor(flag3, levels = c("FM", "UM", "IM"), ordered = TRUE), 
                           alpha = factor(flag3, levels = c("FM", "UM", "IM"), ordered = TRUE)), 
                       bins=100,  position="identity", col = NA) +
        ylab("counts") + ggtitle(sample_id) + 
        geom_text(data = 
                  data.frame(x = 0.75, y = 1E6, tex=paste0("cor=", 
                  dt[flag3 != "NA", round2(cor(m_b, m_f), digits = 3)])), 
                  aes(x=x, y=y, label=tex)) +
        theme(panel.grid.major.x = element_line( size=.1, color="black" ) ) +
        scale_x_continuous(name = "mn - mf,n", breaks=seq(-1, 1, 0.25)) +
        scale_y_continuous(name = "counts", breaks=seq(0, 15E3, 5E3), labels=seq(0, 15, 5)) +
        coord_cartesian(ylim = c(0, 16000)) +
        scale_fill_manual(name="state", values=c("UM" = "blue", "IM" = "darkturquoise", "FM" = "gold2")) +
        scale_alpha_manual(name="state", values=c("UM" = 0.35, "IM" = 0.75, "FM" = 0.35))
      ggsave(plot = p2_zoom, filename = paste0(f_prefix, "_vs_normal_hist_diff_all_zoom.pdf"), device = "pdf",
             width=4, height = 3, unit="in")
      
       p3 <-  ggplot(dt[flag=="UM"|flag2=="UM",], aes(x=m_b-m_f, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1)) 
      ggsave(p3, filename = paste0(f_prefix, "_vs_normal_hist_diff_UM.pdf"), device=pdf,
             width=4, height = 3, unit="in")
      
      p4 <-  ggplot(dt[flag=="FM"|flag2=="FM",], aes(x=m_b-m_f, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1))  
      ggsave(p4, filename = paste0(f_prefix, "_vs_normal_hist_diff_FM.pdf"), device=pdf,
             width=4, height = 3, unit="in")
      
      p5 <-  ggplot(dt[flag=="IM"|flag2=="IM",], aes(x=m_b-m_f, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1)) 
      ggsave(p5, filename = paste0(f_prefix, "_vs_normal_hist_diff_IM.pdf"), device=pdf,
             width=4, height = 3, unit="in")
      
      p6 <-  ggplot(dt[!is.na(DMP_b_corr),], aes(x=m_b-m_f, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1)) 
      ggsave(p5, filename = paste0(f_prefix, "_vs_normal_hist_diff_DMPs.pdf"), device=pdf,
             width=4, height = 3, unit="in")
      
      # annotate clin info
      dt[, sample := sample_id]
      dt[, patient := patient_id]
      dt[, region := region_id]
     
      # save results 
      tmp <- dt[flag!="NA", 
                .(m=mean(m_b-m_f),
                  m_q025=quantile(m_b-m_f, 0.025),
                  m_q975=quantile(m_b-m_f, 0.975),
                  m_q05=quantile(m_b-m_f, 0.05),
                  m_q95=quantile(m_b-m_f, 0.95),
                  m_abs=mean(abs(m_b-m_f)),
                  m_abs_q95=quantile(abs(m_b-m_f), 0.95),
                  m_abs_q99=quantile(abs(m_b-m_f), 0.99),
                  mad=mad(abs(m_b-m_f))), 
                by=.(patient, region, sample, flag)]
      tmp2 <- dt[,
                 .(m=mean(m_b-m_f),
                   m_q025=quantile(m_b-m_f, 0.025),
                   m_q975=quantile(m_b-m_f, 0.975),
                   m_q05=quantile(m_b-m_f, 0.05),
                   m_q95=quantile(m_b-m_f, 0.95),
                   m_abs=mean(abs(m_b-m_f)),
                   m_abs_q95=quantile(abs(m_b-m_f), 0.95),
                   m_abs_q99=quantile(abs(m_b-m_f), 0.99),
                   mad=mad(abs(m_b-m_f))), 
                 by=.(patient, region, sample)]
      tmp2[, flag := "all"]
      dt_list[[length(dt_list)+1]] <- rbind(tmp, tmp2)
      dt_list2[[length(dt_list2)+1]] <- 
        dt[, .(num=length(m_b),
               m=mean(m_b-m_f),
               m_q025=quantile(m_b-m_f, 0.025),
               m_q975=quantile(m_b-m_f, 0.975),
               m_q05=quantile(m_b-m_f, 0.05),
               m_q95=quantile(m_b-m_f, 0.95),
               m_abs=mean(abs(m_b-m_f)),
               m_abs_q95=quantile(abs(m_b-m_f), 0.95),
               m_abs_q99=quantile(abs(m_b-m_f), 0.99),
               mad=mad(abs(m_b-m_f))), by=.(patient, region, sample, flag, flag2)]
      rm(dt, tmp, tmp2, path_output, f_name)
      rm(list=ls(pattern="^p[0-1]"))
      }
}

# rbind all CpGs classed as IM with 
dt2 <- rbindlist(dt_list2) ; rm(dt_list2)
dt2[, flag3 := flag==flag2]
tmp <- dt2[flag != "NA" & flag2 != "NA", .(num=round2(median(num), digits=0)), .(flag,flag2)]
x <- c("UM", "IM", "FM")  
y <- tmp[match(paste(x, rep(x,each=3)),paste(tmp$flag, tmp$flag2)), (num)]
y[is.na(y)] <- 0
mat <- matrix(y, 3,3)
rownames(mat) <- x
colnames(mat) <- x

# tmp <- dt2[flag!="NA"&flag2!="NA", .(num=round2(median(num), digits=0)), .(flag,flag3)]
# x <- c("UM", "IM", "FM")  
# x2 <- c("FALSE", "TRUE") 
# y <- tmp[match(paste(rep(x,2), rep(x2, each=3)),paste(tmp$flag, tmp$flag3)), (num)]
# mat <- matrix(y, 3,2)
# rownames(mat) <- x
# colnames(mat) <- x2
suppressWarnings(rm(x, y, x2))

# rbind all CpGs with highest delta m 
dt3 <- rbindlist(dt_list3) ; rm(dt_list3)
# tmp <- dt3[flag=="IM"|flag2=="IM", 
#            .(l=length(m_b),diff=mean(abs(m_b-m_f))),
#            by=.(chromosome, start, end)]
# setkeyv(tmp, c("chromosome", "start", "end"))
# ggplot(tmp,aes(x=factor(l, 1:5, ordered = TRUE), y=diff))+
#   theme_minimal()+
#   geom_violin(alpha=0.25)

# rbind all CpGs classed as IM with 
dt <- rbindlist(dt_list)
rm(dt_list)

# annotate clin info
tmp <- patient_info[Pop != 1, .(Patient, Region, Pop)]
tmp[, ltx_id := paste0(Patient, "_", Region)]
tmp[, rho := patient_clin[match(tmp$ltx_id, LTX_ID), rho_reseg]]
tmp[, psi := patient_clin[match(tmp$ltx_id, LTX_ID), psi_reseg]]
tmp2 <- dt[, .(x=max(abs(m_q025), abs(m_q975))), by=.(sample, patient, region)]
tmp[, diff := tmp2[match(tmp$Patient, tmp2$patient),.(x)]]
tmp[, effect_max_diff := diff*(1-rho)*2/((1-rho)*2+rho*psi)]

pth_name <- file.path(work_dir, "outputs")
ggplot(dt, aes(x=factor(flag, levels=c("UM", "IM", "FM", "all"), ordered=TRUE), y=m)) + 
  theme_classic() + ylim(c(-0.05, 0.15)) + xlab("Bulk normal methylation state") +
  geom_boxplot() + geom_jitter(width=0.1) + ylab("Bulk - FACS normal methylation") +
  ggsave(filename = paste0(pth_name, "m_diff_normals_by_methylation_state.pdf"), device=pdf,
         width=4, height = 3, unit="in")
ggplot(dt, aes(x=factor(flag, levels=c("UM", "IM", "FM", "all"), ordered=TRUE), y=m_abs)) + 
  theme_classic() + ylim(c(0, 0.22)) + xlab("Bulk normal methylation state") +
  geom_boxplot() + geom_jitter(width=0.1) + ylab("Abs bulk - FACS normal methylation") +
  ggsave(filename = paste0(pth_name, "abs_m_diff_normals_by_methylation_state.pdf"), device=pdf,
         width=4, height = 3, unit="in")
ggplot(dt, aes(x=factor(flag, levels=c("UM", "IM", "FM", "all"), ordered=TRUE), y=mad)) + 
  theme_classic() + ylim(c(0, 0.16)) + xlab("Bulk normal methylation state") +
  geom_boxplot() + geom_jitter(width=0.1) + ylab("Bulk vs. FACS normal methylation dev") +
  ggsave(filename = paste0(pth_name, "deviation_normals_by_methylation_state.pdf"), device=pdf,
         width=4, height = 3, unit="in")

################################################################################
# Plot FACS vs. CAMDAC bulk and purified tumour methylation rate
################################################################################

# set params
per_sample_plots <- FALSE

# loop over each tumour sample
dt_list <- list()
for(sample_id in patient_info[, sort(Sample_Name)]){
  patient_id <- patient_info[Sample_Name == sample_id, Patient]
  region_id <- patient_info[Sample_Name == sample_id, Region]
  pop_id <- paste0("pop", patient_info[Sample_Name == sample_id, Pop])

  if (pop_id == "pop2") {
    print(sample_id)

    # create dirs
    path <- paste0(parent_dir, patient_id, "/")
    path2 <- paste0(parent_dir2, patient_id, "/")

    # get file name prefix
    path_output <- paste0(path, "Methylation/", sample_id, "/")
    path_output2 <- paste0(path2, "Methylation/", region_id, "/re_segmented/")

    # set filenames and filepaths
    output_file <- paste0(path_output, "purified_tumour_DMP_calls.RData")
    output_file2 <- paste0(path_output2, "purified_tumour_DMP_calls.RData")

    # load file processed methylation info
    load(output_file)
    dt1 <- data.table(df_mt)[, .SD, .SDcols = c("chromosome", "start", "end",
                                             "m_b", "m_b_low", "m_b_high",
                                             "DMP_b", "DMP_b_corr")]
    setkeyv(dt1, c("chromosome", "start", "end"))
    colnames(dt1) <- gsub("m_b", "m_f_t", colnames(dt1))
    
    load(output_file2)
    dt2 <- data.table(df_mt)[, .SD, .SDcols=c("chromosome", "start", "end",
                                             "m_n", "m_n_low", "m_n_high",
                                             "m_b", "m_b_low", "m_b_high", 
                                             "m_t_corr", "m_t_low", "m_t_high", 
                                             "DMP_t", "DMP_t_corr")]
    setkeyv(dt2, c("chromosome", "start", "end"))
    
    # merge data.table
    dt <- dt1[dt2, nomatch = 0]

    # subset to DMPs
    dt <- dt[!is.na(DMP_b_corr)]

    # workspace clean up
    suppressWarnings(rm(df_mt, output_file2, output_file, path2, path, dt1, dt2))
    
    # class CpGs
    dt[, flag := ifelse(m_f_t_high < 0.2, "UM", 
                        ifelse(m_f_t_low > 0.8, "FM", 
                               ifelse(m_f_t_low > 0.2 & m_f_t_high < 0.8, "IM", 
                                      "NA")))]
    
    if (per_sample_plots == TRUE) {
      # set the image output file prefix
      f_prefix <- paste0(work_dir, "outputs/", sample_id, "FACS_tumour_vs_CAMDAC_pure_tumour_")

      # plot the methylation histogram
      p1 <- ggplot(dt, aes(y=after_stat(count))) +  theme_classic() + 
        geom_histogram(aes(col = "CAMDAC", fill = "CAMDAC", x = m_t_corr), alpha = 0.1, bins = 25) + 
        geom_histogram(aes(col = "FACS", fill = "FACS", x = m_f_t), alpha = 0.1, bins = 25) + 
        ylab("counts") + ggtitle(sample_id) + xlab("purified tumour methylation rate") +
        scale_color_manual(name = "platform", values = c("CAMDAC" = "orange3", "FACS" = "purple3")) +
        scale_fill_manual(name = "platform", values = c("CAMDAC" = "orange3", "FACS" = "purple3")) + 
        geom_text(data = data.frame(x = 0.15, y = 1E6, 
                                    tex = paste0("cor=", dt[, round2(cor(m_f_t, m_t_corr), digits = 3)])), 
                  aes(x = x, y = y, label = tex)) 
      ggsave(p1, filename = paste0(f_prefix, "_hist_all.pdf"), device = pdf,
            width = 4, height = 3, unit = "in")
      
      # get the quantiles
      q1 <- dt[,quantile(m_t_corr - m_f_t, 0.025)]
      q2 <- dt[,quantile(m_t_corr - m_f_t, 0.975)]
      p2 <- ggplot(dt, aes(x = m_t_corr - m_f_t, y = after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha = 0.1, bins = 100, col = "forestgreen", fill = "forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1)) +
        geom_vline(xintercept = c(q1, q2), linetype="dashed")
      ggsave(p2, filename = paste0(f_prefix, "_hist_diff_all.pdf"), device=pdf,
            width=4, height = 3, unit="in")
      
      p3 <-  ggplot(dt[flag == "UM",], aes(x = m_t_corr - m_f_t, y = after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1)) 
      ggsave(p3, filename = paste0(f_prefix, "_hist_diff_UM.pdf"), device=pdf,
            width=4, height = 3, unit="in")
      
      p4 <-  ggplot(dt[flag == "FM",], aes(x=m_t_corr-m_f_t, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1))
      ggsave(p4, filename = paste0(f_prefix, "_hist_diff_FM.pdf"), device=pdf,
            width=4, height = 3, unit="in")
      
      p5 <-  ggplot(dt[flag == "IM",], aes(x=m_t_corr-m_f_t, y=after_stat(count))) +  theme_classic() + 
        geom_histogram(alpha=0.1, bins=100, col="forestgreen", fill="forestgreen") +
        ylab("counts") + ggtitle(sample_id) + xlim(c(-1, 1))
      ggsave(p5, filename = paste0(f_prefix, "_hist_diff_IM.pdf"), device=pdf,
            width=4, height = 3, unit="in")
    }
    
    dt[, sample := sample_id]
    tmp <- dt[flag != "NA", .(m_t_vs_m_f_t = mean(m_t_corr - m_f_t),
                  m_b_vs_m_f_t = mean(m_b - m_f_t),
                  m_t_vs_m_f_t_abs = mean(abs(m_t_corr - m_f_t)),
                  m_b_vs_m_f_t_abs = mean(abs(m_b - m_f_t)),
                  m_t_vs_m_f_t_mad = mad(abs(m_t_corr - m_f_t)),
                  m_b_vs_m_f_t_mad = mad(abs(m_b - m_f_t))), 
              by=.(sample, flag)]
    tmp2 <- dt[, .(m_t_vs_m_f_t = mean(m_t_corr - m_f_t),
                   m_b_vs_m_f_t = mean(m_b - m_f_t),
                   m_t_vs_m_f_t_abs = mean(abs(m_t_corr - m_f_t)),
                   m_b_vs_m_f_t_abs = mean(abs(m_b - m_f_t)),
                   m_t_vs_m_f_t_mad = mad(abs(m_t_corr - m_f_t)),
                   m_b_vs_m_f_t_mad = mad(abs(m_b - m_f_t))),  
               by=.(sample)]
    tmp2[, flag := "all"]
    dt_list[[length(dt_list) + 1]] <-  rbind(tmp, tmp2)
    rm(dt, tmp, prefix, output_file_DMPs, path_output, f_name)
    rm(list=ls(pattern="^p[0-1]"))
  }
}
dt <- rbindlist(dt_list)
rm(dt_list)

# plot summary metrics across FACS purified tumour samples and samples from the same patient
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_t_vs_m_f_t)) + 
  theme_classic() + ylim(c(-0.3, 0.3)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width=0.1) + ylab("CAMDAC mt - FACS mt")
ggsave(filename = file.path(work_dir, "outputs",
          "m_diff_pure_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_b_vs_m_f_t)) + 
  theme_classic() + ylim(c(-0.3, 0.3)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width=0.1) + ylab("CAMDAC mb - FACS mt")
ggsave(filename = file.path(work_dir, "outputs",
          "m_diff_bulk_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_t_vs_m_f_t_abs)) + 
  theme_classic() + ylim(c(0, 0.35)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width=0.1)  + ylab("|CAMDAC mt - FACS mt|")
ggsave(filename = file.path(work_dir, "outputs",
          "abs_m_diff_pure_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_b_vs_m_f_t_abs)) + 
  theme_classic() + ylim(c(0, 0.35)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width=0.1)  + ylab("|CAMDAC mb - FACS mt|")
ggsave(filename = file.path(work_dir, "outputs",
          "abs_m_diff_bulk_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_t_vs_m_f_t_mad)) + 
  theme_classic() + ylim(c(0, 0.25)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width = 0.1) + ylab("CAMDAC - FACS mt MAD")
ggsave(filename = file.path(work_dir, "outputs",
          "deviation_pure_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")
ggplot(dt, aes(x = factor(flag, levels=c("UM", "IM", "FM", "all"), ordered = TRUE),
               y = m_b_vs_m_f_t_mad)) + 
  theme_classic() + ylim(c(0, 0.25)) + xlab("Tumour methylation state") +
  geom_boxplot() + geom_jitter(width = 0.1)  + ylab("CAMDAC - FACS mt MAD")
ggsave(filename = file.path(work_dir, "outputs",
          "deviation_bulk_tumour_by_methylation_state.pdf"),
       device = pdf, width = 4, height = 3, unit = "in")

################################################################################
# Plot FACS vs. CAMDAC bulk and purified tumour methylation rate
################################################################################

for (sample_id in patient_info[, sort(Sample_Name)]) {
  patient_id <- patient_info[Sample_Name == sample_id, Patient]
  region_id <- patient_info[Sample_Name == sample_id, Region]
  pop_id <- paste0("pop", patient_info[Sample_Name == sample_id, Pop])

  if (pop_id == "pop2") {
    print(sample_id)
    
    # create dirs
    path <- paste0(parent_dir, patient_id, "/")
    path2 <- paste0(parent_dir2, patient_id, "/")
    
    # get file name prefix
    path_output <-paste0(path, "Methylation/", sample_id, "/")
    path_output2 <-paste0(path2, "Methylation/", region_id, "/re_segmented/")
    
    # set filenames and filepaths
    output_file = paste0(path_output, "purified_tumour_DMP_calls.RData")
    output_file2 = paste0(path_output2, "purified_tumour_DMP_calls.RData")
    
    # load FACS processed methylation info
    load(output_file)
    dt1 <- data.table(df_mt)[, .SD, .SDcols=c("chromosome", "start", "end", 
                                              "m_b", "m_b_low", "m_b_high", 
                                              "m_n", "m_n_low", "m_n_high",
                                              "DMP_b", "DMP_b_corr")]
    setkeyv(dt1, c("chromosome", "start", "end"))

    # Replace both 'm_b' and 'm_n' patterns in column names
    colnames(dt1) <- gsub("m_b", "m_t", colnames(dt1))
    colnames(dt1) <- gsub("m_(t|n)", "m_f_\\1", colnames(dt1))

    # load FACS processed methylation info
    load(output_file2)
    dt2 <- data.table(df_mt)[, .SD, .SDcols=c("chromosome", "start", "end", 
                                              "m_n", "m_n_low", "m_n_high",
                                              "m_b", "m_b_low", "m_b_high",
                                              "m_t_corr", "m_t_low", "m_t_high", 
                                              "DMP_t", "DMP_t_corr")]
    setkeyv(dt2, c("chromosome", "start", "end"))
    
    # merge data.table (remove rows without info in both tables)
    dt <- dt1[dt2, nomatch = 0]
    suppressWarnings(rm(df_mt, output_file2, output_file, path2, path, dt1, dt2))
    
    # plot 2d density (FACS purified vs. CAMDAC purified tumour)
    df_cor <- data.frame(
      x = 0.15,
      y = 0.9,
      tex = paste0("cor=", dt[!is.na(m_f_t), round2(cor(m_f_t, m_t_corr), digits=3)])
    )
    p <- ggplot(dt, aes(x = m_f_t, y = m_t_corr)) +  theme_classic() + 
        geom_bin2d(bins=50) + 
        scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 2, breaks = c(10, 100, 1E3, 1E4, 1E5), 
                             labels = parse(text=c("10^1","10^2","10^3","10^4","10^5")),
                             na.value = "white", limits = c(10, 10^(4.5)), oob = squish) +
        geom_text(data = df_cor, aes(x = x, y = y, label = tex)) +
        guides(colorbar=guide_legend(title=expression(count~(log[10])))) +
        ylab("m_t CAMDAC") + ggtitle(sample_id) + xlab(paste0("m_t FACS"))
      
    ggsave(p, filename = paste0(work_dir, "/outputs/FACS_", sample_id, "_vs_pure_tumour.pdf"), 
          device = pdf, width = 4, height = 3, unit = "in")
    rm(p, df_cor)

    # plot 2d density (FACS purified vs. bulk tumour)
    df_cor <- data.frame(
      x = 0.15,
      y = 0.9,
      tex = paste0("cor=", dt[!is.na(m_f_t), round2(cor(m_f_t, m_b), digits = 3)])
      )
    p <- ggplot(dt, aes(x = m_f_t, y = m_b)) +  theme_classic() + 
        geom_bin2d(bins=50) + 
        scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 2, breaks = c(10, 100, 1E3, 1E4, 1E5), 
                             labels = parse(text=c("10^1","10^2","10^3","10^4","10^5")),
                             na.value = "white", limits = c(10, 10^(4.5)), oob = squish) +
        geom_text(data = df_cor, aes(x = x, y = y, label = tex)) +
        guides(colorbar=guide_legend(title=expression(count~(log[10])))) +
        ylab("m_b") + ggtitle(sample_id) + xlab(paste0("m_t FACS"))
      
    ggsave(p, filename = paste0(work_dir, "/outputs/FACS_", sample_id, "_vs_bulk_tumour.pdf"), 
          device = pdf, width = 4, height = 3, unit = "in")
    rm(p, df_cor)
      
    # plot 2d density (FACS purified vs. CAMDAC purified tumour methylation rates at DMPs)
    df_cor <- data.frame(
      x = 0.15,
      y = 0.9,
      tex = paste0("cor=",
        dt[!is.na(m_f_t) & !is.na(DMP_b_corr), round2(cor(m_f_t, m_t_corr), digits = 3)]
      )
    )
    p <- ggplot(dt[!is.na(DMP_b_corr)], aes(x = m_f_t, y = m_t_corr)) +  theme_classic() + 
        geom_bin2d(bins=50) + 
        scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 1.1, breaks = c(1, 10, 100, 1E3), 
                             labels = parse(text=c("10^0", "10^1","10^2","10^3")),
                             na.value = "white", limits = c(1, 10^(2.5)), oob = squish) +
        geom_text(data = df_cor, aes(x = x, y = y, label = tex)) +
        guides(colorbar=guide_legend(title=expression(count~(log[10])))) +
        ylab("m_t CAMDAC") + ggtitle(sample_id) + xlab(paste0("m_t FACS"))
      
    ggsave(p, filename = paste0(work_dir, "/outputs/FACS_", sample_id, "_vs_pure_tumour_DMP.pdf"), 
          device = pdf, width = 4, height = 3, unit = "in")
    rm(p, df_cor)
      
    # plot 2d density (FACS purified vs. CAMDAC bulk tumour methylation rates at DMPs)
    df_cor <- data.frame(
      x = 0.15,
      y = 0.9,
      tex = paste0("cor=",
        dt[!is.na(m_f_t) & !is.na(DMP_b_corr), round2(cor(m_f_t, m_b), digits = 3)]
      )
    )
    p <- ggplot(dt[!is.na(DMP_b_corr)], aes(x = m_f_t, y = m_b)) +  theme_classic() + 
        geom_bin2d(bins=50) + 
        scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 1.1, breaks = c(1, 10, 100, 1E3), 
                             labels = parse(text=c("10^0", "10^1","10^2","10^3")),
                             na.value = "white", limits = c(1, 10^(2.5)), oob = squish) +
        geom_text(data = df_cor, aes(x = x, y = y, label = tex)) +
        guides(colorbar=guide_legend(title=expression(count~(log[10])))) +
        ylab("m_b") + ggtitle(sample_id) + xlab(paste0("m_t FACS"))
      
    ggsave(p, filename = paste0(work_dir, "/outputs/FACS_", sample_id, "_vs_bulk_tumour_DMP.pdf"), 
          device = pdf, width = 4, height = 3, unit = "in")
    rm(p, df_cor)
  }
} 
