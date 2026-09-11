#' @description To showcase CAMAC DMR calling, we evaluate recall at
#' well-known LUAD and LUSC promoter DMRs and compare with bulk calls.
#' 
#' @author Lili Larose Cadieux

# Load library
x <- c("stringr", "GenomicRanges", "data.table", "fst", "ggplot2",
       "gridExtra", "svglite")

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
# Set paths and load inputs
################################################################################

# set parent directory
tx_dir <- "~/rdeocd/txMeth"

# set the working directory
work_dir <- "~/data/CAMDAC_paper/figure3FG"
setwd(work_dir)

# set the output directory
if(!dir.exists("outputs")) dir.create("outputs")

# load patient_info
f_nm <-
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
patient_info <- data.table::fread(file.path(tx_dir, f_nm))
rm(f_nm)

# remove lymph nodes
patient_info <- patient_info[!grepl("LN", Sample)]
patient_ids <- unique(patient_info$Patient)

# set file prefixes - we will be focusing on promoter DMRs
preprefix <- "_promoter"
prefixs <- c("bulk_tumour", "purified_tumour")

# load the genes of interest
belinksky_genes <- data.table::fread("inputs/belinsky_gene_promoters.tsv")

# load the prevalence
belinksky_genes_prevalence <- data.table::fread("inputs/belinsky_gene_prevalence.tsv")

# excluded genes:
# CASP8 - not found to be hyper in NSCLC
# FAS and TRAIL1 are only hyper in SCLC
# refs - see thesis page 103:
# + GATA4, GATA5, PAX5α and PAX5β # 10.1158/1940-6207.CAPR-17-0177
# 173 - Belinsky, Nature Reviews Cancer, 4(9):707–717, 2004
# 169 - Tang, Journal of the National Cancer Institute, 92(18):1511–1516, 2000
# 171 - Kim, Cancer research, 61(8):3419–24, 2001
# 175 - Kim, Oncogene, 20(14):1765–1770, 2001
# 176 - Zöchbauer-Müller, Cancer Research, 2001
# 177 - Bearzatto, Clinical cancer research, 8(12):3782–7, 2002.
# 178 - Palmisano, Cancer research, 63(15):4620–5, 2003.
# 179 - Mizuno, Oncogene, 21(15):2328–2333, 2002.
# 180 - Corn, Carcinogenesis, 24(1):47–51, 2003.
# 181 - Tomizawa, Lung Cancer, 46(3):305–312, 2004.
# 170 - Seng, British Journal of Cancer, 99(2):375–382, 2008.
# 182 - Esteller, Cancer research, 59(4):793–7, 1999.
# 183 - Toyooka, Cancer research, 61(11):4556–60, 2001.
# 184 - Sathyanarayana, Clinical cancer research, 9(7):2665–2672, 2003.
# 185 - Dammann, Nature Genetics, 25(3):315–319, 2000.
# 186 - Burbee, Journal of the National Cancer Institute, 93(9):691–9, 2001.
# 187 - Hesson, Oncogene, 22(6):947–954, 2003.

################################################################################
# Gather promoters with DMPs
################################################################################

# Look for DMPs at recurrently DM gene promoters from lit. review
belinksky_genes$camdac_region_id <- character()
new_cols <-
  c("num_CpGs_at_least_2_samples", "num_CpGs_avg_per_sample",
    "num_samples_at_least_3_CpGs",
    "cov_t_avg_per_sample_inc_0", "cov_t_avg_per_sample_exc_0",
    "num_DMPs_avg_per_sample", "num_DMPs_max_per_sample",
    "num_hyper_avg_per_sample", "num_hyper_max_per_sample")
belinksky_genes[,
  (new_cols) :=
    replicate(9, numeric(), simplify = FALSE)
  ]

# load CAMDAC RRBS annotations
f_nm <- "hg19_all_regions_annotations.fst"
anno_dir <- "~/rdeocd/pipeline_files/hg19_annotations"
camdac_regions <-
  fst::read_fst(path = file.path(anno_dir, f_nm), as.data.table = TRUE)
rm(anno_dir, f_nm)

# set columns of interest from the CAMDAC methylation outputs
id_cols <- c("patient", "sample", "ID", "chromosome", "start", "end")
prom_cols <- c("promoter", "CGI")
cov_cols <- paste0("cov_", c("n", "b", "t"))
m_cols <- paste0("m_", c("n", "b", "t"))
#diff_cols <- paste0("m_diff_", c("tn", "bn"))
dmp_cols <- paste0("DMP_", c("tn", "bn"))
all_cols <- c(id_cols, cov_cols, m_cols, dmp_cols, prom_cols)

# load all CpGs covered in at least 2 samples
all_cpgs_f_nm <- file.path(tx_dir, "CpG_wise_CAMDAC_results_all.fst")
dt_cpgs <-
  fst::read_fst(path = all_cpgs_f_nm, columns = all_cols, as.data.table = TRUE)

gene_methylation_lists <- list()
for (i in seq_len(nrow(belinksky_genes))) {
  # get the gene and promoter ID
  gene_name <- belinksky_genes$gene[i]
  gene_prom_id <- paste0(gene_name, belinksky_genes$promoter_id[i])

  # get the coordinates
  gene_chrom <- belinksky_genes$chrom[i]
  gene_start <- belinksky_genes$start[i]
  gene_end <- belinksky_genes$end[i]

  # load ref, is prom is anno?
  dt_regions <-
    camdac_regions[chrom == gene_chrom &
            start - 1 <= gene_end & # -1 for CCGG loci
            end + 1 >= gene_start] # + 1 for CCGG loci
  c_id <- dt_regions$cluster_id
  belinksky_genes$camdac_region_id[i] <-
    paste(c_id, collapse = ",")

  # get the CpGs of interest
  dt_gene <-
    dt_cpgs[chromosome == gene_chrom &
            start - 1 >= gene_start & # -1 for CCGG loci
            end + 1 <= gene_end] # + 1 for CCGG loci

  if (nrow(dt_gene) == 0) {
    belinksky_genes[i, (new_cols) := 0]
  } else {
    # annotate gene id and save gene object
    dt_gene[, gene_id := gene_prom_id]
    gene_methylation_lists[[length(gene_methylation_lists) + 1]] <-
      dt_gene

    # subset to tumour samples
    dt_gene_t <- dt_gene[sample != "N"]

    # how many unique CpGs are covered in any 2 samples
    num_CpGs_total <- dt_gene_t[, length(unique(ID))]
    belinksky_genes[i, (new_cols[1]) := num_CpGs_total]

    # how many unique CpGs are covered per sample
    num_samples_total <- dt_gene_t[, length(unique(paste(patient, sample)))]
    num_cpgs_per_sample <-
      dt_gene_t[, .(num_cpgs = length(ID)), by = .(patient, sample)]$num_cpgs
    
    if (num_samples_total == 122) {
      avg_num_cpgs_per_sample <- mean(num_cpgs_per_sample)
    } else {
      samples_without_cov <- base::rep(0, times = 122 - num_samples_total)
      avg_num_cpgs_per_sample <- mean(c(samples_without_cov, num_cpgs_per_sample))
    }
    belinksky_genes[i, (new_cols[2]) := avg_num_cpgs_per_sample]

    # get the number of samples with at least 3 CpGs covered (needed to have DMR)
    belinksky_genes[i, (new_cols[3]) := sum(num_cpgs_per_sample >= 3, na.rm = TRUE)]

    # note down the tumour coverage
    cov_t_per_sample <-
      dt_gene_t[, .(cov_t = mean(cov_t, na.rm = TRUE)), by = .(patient, sample)]$cov_t

    if (num_samples_total == 122) {
      cov_t_per_sample_inc_0 <- mean(cov_t_per_sample)
    } else {
      cov_t_per_sample_inc_0 <- mean(c(samples_without_cov, cov_t_per_sample))
    }
    belinksky_genes[i, (new_cols[4]) := cov_t_per_sample_inc_0]
    belinksky_genes[i, (new_cols[5]) := mean(cov_t_per_sample, na.rm = TRUE)]

    # how many unique DMPs per sample
    num_dmps_per_sample <-
      dt_gene_t[, .(num_dmps = sum(!is.na(DMP_tn))), by = .(patient, sample)]$num_dmps
    belinksky_genes[i, (new_cols[7]) := max(num_dmps_per_sample)]

    if (num_samples_total == 122) {
      avg_num_dmps_per_sample <- mean(num_dmps_per_sample)
    } else {
      avg_num_dmps_per_sample <- mean(c(samples_without_cov, num_dmps_per_sample))
    }
    belinksky_genes[i, (new_cols[6]) := avg_num_dmps_per_sample]

    # how many hypermethylated loci per sample
    num_hyper_per_sample <-
      dt_gene_t[,
        .(num_hyper = sum(!is.na(DMP_tn) & DMP_tn == "hyper")),
        by = .(patient, sample)]$num_hyper
    belinksky_genes[i, (new_cols[7]) := max(num_hyper_per_sample)]

    if (num_samples_total == 122) {
      avg_num_hyper_per_sample <- mean(num_hyper_per_sample)
    } else {
      avg_num_hyper_per_sample <- mean(c(samples_without_cov, num_hyper_per_sample))
    }
    belinksky_genes[i, (new_cols[8]) := avg_num_hyper_per_sample]
  }
}

# save outputs
data.table::fwrite(belinksky_genes, "outputs/dmps_at_belinsky_gene_proms.tsv")
saveRDS(gene_methylation_lists,
  file = "outputs/cpg_methylation_belinsky_genes.rds"
  )

################################################################################
# DMR recall analysis
################################################################################

# get all genes with a EPD promoter covered by at least 3 CpGs
belinksky_genes <- fread("outputs/dmps_at_belinsky_gene_proms.tsv")
belinksky_genes[, table(num_samples_at_least_3_CpGs)]
belinksky_genes[num_samples_at_least_3_CpGs > 1, unique(gene)]
belinksky_genes_subset <- belinksky_genes[num_samples_at_least_3_CpGs > 1]

# get unique patient ids
patients <- patient_info[, unique(Patient)]
patients_luad <- patient_info[grep("adeno", Histology), unique(Patient)]
patients_lusc <-
  patient_info[grep("adeno", Histology, invert = TRUE), unique(Patient)]

# create results object
region_bulk_methylation_list <- list()
region_pure_methylation_list <- list()

# iterate over patients
for (i in seq_len(nrow(patient_info))) {
  # get patient and sample id
  patient <- patient_info[i, Patient]
  sample <- patient_info[i, Sample]

  # iterate over bulk and purified tumuour DMR data
  for (prefx in prefixs) {
    # load patient object
    dir_nm <- file.path(tx_dir, patient, "Methylation", sample, "re_segmented")
    f_nm <- paste0(prefx, "_promoter_DMR_calls.RData")
    if (!file.exists(file.path(dir_nm, f_nm))) {
      dir_nm <- file.path(tx_dir, patient, "Methylation", sample, "raw")
    }
    if (!file.exists(file.path(dir_nm, f_nm))) {
      dir_nm <- file.path(tx_dir, patient, "Methylation", sample)
    }
    load(file.path(dir_nm, f_nm))
    dt_dmrs <- df_annotated_regions

    # iterate over genes
    for (j in seq_len(nrow(belinksky_genes_subset))) {
      # subset
      dt_dmrs_gene <-
        dt_dmrs[cluster_id == belinksky_genes_subset$camdac_region_id[j]]

      # annotate patient and sample ID
      dt_dmrs_gene$patient <- patient
      dt_dmrs_gene$sample <- sample

      # annotate gene
      dt_dmrs_gene$gene <- belinksky_genes_subset$gene[j]
      dt_dmrs_gene$promoter_id <- belinksky_genes_subset$promoter_id[j]

      # store data
      if (prefx == "bulk_tumour") {
        region_bulk_methylation_list[[
          length(region_bulk_methylation_list) + 1
          ]] <- dt_dmrs_gene
      } else {
        region_pure_methylation_list[[
          length(region_pure_methylation_list) + 1
          ]] <- dt_dmrs_gene
      }
    }
  }
}

# combine methylation data
region_bulk_methylation <- rbindlist(region_bulk_methylation_list)
region_pure_methylation <- rbindlist(region_pure_methylation_list)
region_bulk_methylation$method <- "bulk_tumour"
region_pure_methylation$method <- "purified_tumour"

# compute the methylation difference
region_pure_methylation[, m_diff_DMPs :=
  as.numeric(m_t_DMPs) - as.numeric(m_n_DMPs)]
region_bulk_methylation[, m_diff_DMPs :=
  as.numeric(m_b_DMPs) - as.numeric(m_n_DMPs)]

# columns to select
keep_cols <-
  c("patient", "sample", "gene", "promoter_id", "method", "cov_n",
    "m_n", "m_n_DMPs", "m_diff_DMPs", "DMR", "DMP_counts", "consec_DMPs")

# merge data
dt_regions_merge <-
  rbind(
    region_pure_methylation[, .SD, .SDcols = c(keep_cols)],
    region_bulk_methylation[, .SD, .SDcols = c(keep_cols)]
    )

# format histology, LTX and CRUKIDs
dt_regions_merge[, hist := patient_info[match(dt_regions_merge$patient, Patient), Histology]]
dt_regions_merge[, CRUKID := patient_info[match(dt_regions_merge$patient, Patient), CRUKID]]
dt_regions_merge[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
dt_regions_merge[, LTX_ID := paste0(patient, "_", sample)]
dt_regions_merge[, CRUK_IDs := paste0(CRUKID, "_", sample)]

# note down DMR type
dt_regions_merge[, DMR_type :=
  ifelse(sum(is.na(m_diff_DMPs)) == length(m_n), as.character(NA),
      ifelse(sum(m_diff_DMPs > 0, na.rm = TRUE) > 0,
      "hyper", "hypo")),
    by = .(patient, sample, gene, promoter_id, method)]

# summarise across multiple promoters of the same genes
dt_summary_by_gene <-
  dt_regions_merge[, .(
    m_diff_DMPs = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMR = ifelse(sum(!is.na(DMR)) > 0, "DMR", NA_character_),
    DMR_type =
      ifelse(sum(!is.na(DMR) & DMR_type == "hypo") > 0 &
        mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE) < 0,
        "hypo",
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper") > 0 &
        mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE) > 0,
        "hyper", 
      NA_character_))
    ),
  by = .(gene, method, hist, patient, sample)]

# summarise DMRs by method across promoters of the same genes and samples.
dt_summary <-
  dt_summary_by_gene[, .(
    mean_m_diff = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),

    # NSCLC per sample
    DMRs_ps = sum(!is.na(DMR)),
    hypo_ps = sum(!is.na(DMR) & DMR_type == "hypo"),
    hyper_ps = sum(!is.na(DMR) & DMR_type == "hyper")
    ),
  by = .(gene, method, hist)]

# get the differences between bulk and purified data
histologies <- c("LUAD", "LUSC")
number_of_samples <- c(68, 54)
names(number_of_samples) <- histologies
for (hs in histologies) {
  dt_plot <-
    dt_summary[hist == hs,
      .(total = hyper_ps[method == "bulk_tumour"],
        diff = hyper_ps[method == "bulk_tumour"]
              - hyper_ps[method == "purified_tumour"]),
      by = .(gene, hist)]
  
  if (hs == "LUAD") {
    gene_ord <- dt_plot[order(diff, total), gene]
  }
  dt_plot[, gene_name := factor(gene, levels = gene_ord, ordered = TRUE)]
  dt_plot[, ymin := ifelse(diff > 0, 0, diff)]
  dt_plot[, ymax := ifelse(diff < 0, 0, diff)]

  # make plot
  p <- ggplot(dt_plot) +
      geom_linerange(aes(x = gene_name, ymin = 0, ymax = total)) +
      geom_linerange(aes(x = gene_name, ymin = ymin, ymax = ymax), col = "firebrick") +
      geom_point(aes(x = gene_name, total)) +
      geom_point(aes(x = gene_name, diff), col = "firebrick") +
      geom_hline(yintercept = 0) +
      scale_y_continuous("samples") +
      facet_wrap(. ~ hist) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
      
  assign(paste0("p_", tolower(hs)), p)
}

# save plot
f_nm <- file.path("outputs", "DMR_recall_epb.svg")
gd <- arrangeGrob(p_luad, p_lusc, nrow = 1)
ggsave(plot = gd, filename = f_nm, device = "svg", width = 6, height = 2.75)

################################################################################
# DMRs recall notes - EPD naive
################################################################################

# create results object
summary_list <- list()

# columns to select
keep_cols <-
  c("patient", "sample", "external_gene_name", "method", 
    "cov_n", "m_n", "m_n_DMPs", "m_diff_DMPs", "DMR",
    "CpG_counts", "DMP_counts", "consec_DMPs")

# iterate over bulk and purified tumuour DMR data
for (prefx in prefixs) {

  # set DMR file name
  if (prefx == "purified_tumour") {
      f_nm <-
        file.path("outputs", "promoter_regions_camdac_DMR_calls_belinsky.csv")
      # file.path(tx_dir, "promoter_regions_camdac_DMR_calls_belinsky.RData")
  } else {
      f_nm <-
        file.path("outputs", "promoter_regions_bulk_DMR_belinsky.csv")
      # file.path(tx_dir, "promoter_regions_bulk_DMR_belinsky.RData")
  }

  if (! file.exists(f_nm)) {
    # load all genes
    load("inputs/hg19_external_gene_id_UCSC.RData")
    ucsc_gene_names <- genes$external_gene_name
    ucsc_gene_ids <- genes$ensembl_gene_id
    genes <-
      c("CDH1", "CDH13", "CDKN2A", "CHFR", "DAPK1", "DLEC1", "GATA4", "GATA5", "LAMA3",
        "LAMB3", "LAMC3", "MGMT", "MLH1", "PAX5", "RARB", "RASSF1", "RASSF5", "UBE2N",
        "TIMP3")
    ensembl_ids <- ucsc_gene_ids[match(genes, ucsc_gene_names)]
    
    # create results object
    region_methylation_list <- list()

    # iterate over patients
    for (i in seq_len(nrow(patient_info))) {

      # get patient and sample id
      patient <- patient_info[i, Patient]
      sample <- patient_info[i, Sample]
      cat("Processing.... ", patient, " - Sample: ", sample, ".\n", sep = "")

      # load patient promoter region methylation object
      dir_nm <- file.path(tx_dir, patient, "Methylation", sample, "re_segmented")
      patient_f_nm <- paste0(prefx, "_promoter_DMR_calls.RData")
      if (!file.exists(file.path(dir_nm, patient_f_nm))) {
        dir_nm <- file.path(tx_dir, patient, "Methylation", sample, "raw")
      }
      if (!file.exists(file.path(dir_nm, patient_f_nm))) {
        dir_nm <- file.path(tx_dir, patient, "Methylation", sample)
      }
      load(file.path(dir_nm, patient_f_nm))
      dt_tmp_proms <- df_annotated_regions
      rm(df_annotated_regions)

      # get the promoter methylation bins for genes of interest
      dt_proms_list <- list()
      for (j in seq_len(length(genes))) {
        tmp <- dt_tmp_proms[grepl(genes[j], gene_prom_ids), ]
        if (nrow(tmp) > 0) {
          expr <- paste0("(;", genes[j], ",)|(^", genes[j], ",)")
          tmp <- tmp[grepl(expr, gene_prom_ids), ]
        } else {
          tmp <- tmp[1, ]
        }
        tmp$external_gene_name <- genes[j]
        tmp$ensembl_gene_id <- ensembl_ids[j]
        dt_proms_list[[j]] <- tmp
      }
      for (k in length(dt_proms_list):1) {
        if (is.na(dt_proms_list[[k]]$cluster_id[1])) {
          dt_proms_list[[k]] <- NULL
        }
      }

      # combine methylation data across genes of interest
      dt_proms <- rbindlist(dt_proms_list)
      cat("We found ", nrow(dt_proms), " promoter regions in sample ", sample, ".\n", sep = "")
      cat("We found ", nrow(dt_proms[!is.na(DMR)]), " promoter DMRs in ", sample, ".\n", sep = "")

      # annotate patient and sample
      dt_proms$patient <- patient
      dt_proms$sample <- sample
      region_methylation_list[[length(region_methylation_list) + 1]] <-
          dt_proms
      rm(dt_tmp_proms, dt_proms_list, tmp, k, j)
    }

    # save DMRs
    dt_proms_final <- rbindlist(region_methylation_list)
    fwrite(dt_proms_final, file = f_nm, sep = ",")

    # compute the methylation difference
    if (prefx == "purified_tumour") {
      dt_proms_final[, m_diff_DMPs :=
        as.numeric(m_t_DMPs) - as.numeric(m_n_DMPs)]
    } else {
      dt_proms_final[, m_diff_DMPs :=
        as.numeric(m_b_DMPs) - as.numeric(m_n_DMPs)]
    }
    
    # annotate method
    dt_proms_final$method <- prefx

    # add to list
    summary_list[[length(summary_list) + 1]] <-
      dt_proms_final[, .SD, .SDcols = keep_cols]
  } else {
    # load DMRs
    dt_proms <- data.table::fread(file = f_nm)

    # compute the methylation difference
    if (prefx == "purified_tumour") {
      dt_proms[, m_diff_DMPs := as.numeric(m_t_DMPs) - as.numeric(m_n_DMPs)]
    } else {
      dt_proms[, m_diff_DMPs := as.numeric(m_b_DMPs) - as.numeric(m_n_DMPs)]
    }

    # add to list
    dt_proms$method <- prefx
    summary_list[[length(summary_list) + 1]] <-
      dt_proms[, .SD, .SDcols = c(keep_cols)]
  }
}

# bind the results
dt_proms <- rbindlist(summary_list)

# format histology, LTX and CRUKIDs
dt_proms[, hist := patient_info[match(dt_proms$patient, Patient), Histology]]
dt_proms[, CRUKID := patient_info[match(dt_proms$patient, Patient), CRUKID]]
dt_proms[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
dt_proms[, LTX_ID := paste0(patient, "_", sample)]
dt_proms[, CRUK_IDs := paste0(CRUKID, "_", sample)]

# note down DMR type
dt_proms[, DMR_type :=
  ifelse(sum(is.na(DMR)) == length(m_n), as.character(NA),
      ifelse(sum(m_diff_DMPs[!is.na(DMR)] > 0, na.rm = TRUE) > 0,
      "hyper", "hypo")),
    by = .I]

# summarise across multiple genomic bins by genes
dt_summary_by_gene <-
  dt_proms[, .(
    m_diff_DMPs = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMR = ifelse(sum(!is.na(DMR)) > 0, "DMR", NA_character_),
    DMR_type =
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper") > 0,
        "hyper",
      ifelse(sum(!is.na(DMR) & DMR_type == "hypo") > 0,
        "hypo",
      NA_character_)),
    CpG_counts = ifelse(sum(!is.na(DMR)) > 0,
      max(CpG_counts[!is.na(DMR)], na.rm = TRUE),
      max(CpG_counts, na.rm = TRUE)),
    DMP_counts = ifelse(sum(!is.na(DMR)) > 0,
      as.integer(max(DMP_counts[!is.na(DMR)], na.rm = TRUE)),
      NA_integer_),
    consec_DMPs = ifelse(sum(!is.na(DMR)) > 0,
      as.integer(max(consec_DMPs[!is.na(DMR)], na.rm = TRUE)),
      NA_integer_)
    ),
  by = .(external_gene_name, method, patient, sample, hist)]

# summarise DMRs by method across cluster of the same genes and samples.
dt_summary_across_samples <-
  dt_summary_by_gene[, .(
    mean_m_diff = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),

    # NSCLC per sample
    num_samples = length(DMR),
    DMRs_ps = sum(!is.na(DMR)),
    hypo_ps = sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE),
    hyper_ps = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
    ),
  by = .(external_gene_name, method, hist)]

# get the differences between bulk and purified data
histologies <- c("LUAD", "LUSC", "NSCLC")
for (hs in histologies) {
  # subset to histology of interest
  if (hs != "NSCLC") {
    dt_plot <- dt_summary_across_samples[hist == hs]
  } else {
    dt_plot <- dt_summary_across_samples
    dt_plot$hist <- "NSCLC"
  }

  # generate data from dot plot
  dt_plot <-
    dt_plot[,
      .(total_bulk = sum(hyper_ps[method == "bulk_tumour"]),
        total_tumour = sum(hyper_ps[method == "purified_tumour"]),
        diff = sum(hyper_ps[method == "bulk_tumour"]) -
               sum(hyper_ps[method == "purified_tumour"])),
      by = .(external_gene_name, hist)]
  
  # add expected prevalence just to remove genes without reported prevalence
  dt_plot[, prevalence :=
    belinksky_genes_prevalence[histology == hs][
      match(dt_plot$external_gene_name, gene), prevalence]
      ]
  dt_plot <- dt_plot[!is.na(prevalence)]
  
  # order
  gene_ord <- dt_plot[order(total_bulk - diff), external_gene_name]
  dt_plot[, gene_name := factor(external_gene_name, levels = gene_ord, ordered = TRUE)]
  dt_plot[, ymin := ifelse(diff > 0, 0, diff)]
  dt_plot[, ymax := ifelse(diff < 0, 0, diff)]

  # make plot
  p <- ggplot(dt_plot) +
      geom_linerange(aes(x = gene_name, ymin = 0, ymax = total_bulk)) +
      geom_linerange(aes(x = gene_name, ymin = ymin, ymax = ymax), col = "firebrick") +
      geom_point(aes(x = gene_name, total_bulk)) +
      geom_point(aes(x = gene_name, diff), col = "firebrick") +
      geom_hline(yintercept = 0) +
      scale_y_continuous("samples") +
      coord_cartesian(ylim = c(-32, 65)) +
      facet_wrap(. ~ hist) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))

  if (hs != "NSCLC") p <- p + theme(legend.position = "none")
  assign(paste0("p_", tolower(hs)), p)
}

# save plot
gd <- arrangeGrob(p_luad, p_lusc, p_nsclc, nrow = 1, widths = c(1, 1, 0.6))
f_nm <- file.path("outputs", "DMR_recall_camdac_ps")
ggsave(plot = gd, filename = paste0(f_nm, ".svg"), device = "svg", width = 6, height = 2.75)
ggsave(plot = gd, filename = paste0(f_nm, ".pdf"), device = "pdf", width = 6, height = 2.75)

# summarise across multiple samples by patient
dt_summary_by_patient <-
  dt_summary_by_gene[, .(
    m_diff_DMPs = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMR = ifelse(sum(!is.na(DMR)) > 0, "DMR", NA_character_),
    DMR_type =
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) > 0,
        "hyper",
      ifelse(sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE) > 0,
        "hypo",
      NA_character_)),
    CpG_counts = ifelse(sum(!is.na(DMR)) > 0,
      max(CpG_counts[!is.na(DMR)], na.rm = TRUE),
      max(CpG_counts, na.rm = TRUE))
    #DMP_counts = ifelse(sum(!is.na(DMR)) > 0,
    #  as.integer(max(DMP_counts[!is.na(DMR)], na.rm = TRUE)),
    #  NA_integer_),
    #consec_DMPs = ifelse(sum(!is.na(DMR)) > 0,
    #  as.integer(max(consec_DMPs[!is.na(DMR)], na.rm = TRUE)),
    #  NA_integer_)
    ),
  by = .(external_gene_name, method, patient, hist)]

# summarise DMRs by method across cluster of the same genes and patients
dt_summary_across_patient <-
  dt_summary_by_patient[, .(
    mean_m_diff = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    num_patients = length(DMR),
    DMRs_pp = sum(!is.na(DMR)),
    hypo_pp = sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE),
    hyper_pp = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
    ),
by = .(external_gene_name, method, hist)]

# get the differences between bulk and purified data
histologies <- c("LUAD", "LUSC", "NSCLC")
metrics_list <- list()
diff_list <- list()
for (hs in histologies) {
  # subset to histology of interest
  if (hs != "NSCLC") {
    dt_plot <- dt_summary_across_patient[hist == hs]
  } else {
    dt_plot <- dt_summary_across_patient
    dt_plot$hist <- "NSCLC"
  }

  # generate data from dot plot
  dt_plot <-
    dt_plot[,
      .(shared = sum(hyper_pp[method == "bulk_tumour"]) /
              sum(num_patients[method == "purified_tumour"]),
        camdac_only = (sum(hyper_pp[method == "purified_tumour"])
              - sum(hyper_pp[method == "bulk_tumour"])) /
              sum(num_patients[method == "purified_tumour"])),
      by = .(external_gene_name, hist)]

  if (hs == "NSCLC") {
    print(dt_plot[, min(camdac_only)])
    print(dt_plot[, max(camdac_only)])
    print(dt_plot[, mean(camdac_only)])
  }
  dt_melt <- melt(dt_plot, id.vars = c("external_gene_name", "hist"))
  dt_melt[variable == "camdac_only" & value < 0, variable := "bulk_only"]

  # get expected prevalence
  dt_melt[, prevalence :=
    belinksky_genes_prevalence[histology == hs][
      match(dt_melt$external_gene_name, gene), prevalence]
      ]

  # remove genes without reported prevalence
  dt_melt <- dt_melt[!is.na(prevalence)]

  # get range of value where available
  dt_melt[, c("min", "max") := tstrsplit(prevalence, split = "-")]
  dt_melt[, c("min", "max") := lapply(.SD, function(x) as.numeric(x)/100),
          .SDcols = c("min", "max")]
  dt_melt[is.na(max), "max" := min]

  # order genes
  gene_ord <- dt_melt[, .(total = sum(value)), by = .(external_gene_name)]
  gene_ord <- gene_ord[order(total), external_gene_name]
  assign(paste0("gene_ord_", tolower(hs)), gene_ord)
  dt_melt[, gene_name := factor(external_gene_name, levels = gene_ord, ordered = TRUE)]

  # factorise colouring
  # levs <- c("bulk_only", "camdac_only", "shared")
  levs <- c("camdac_only", "shared")
  dt_melt[, variable := factor(variable, levels = levs, ordered = TRUE)]

  # store differences
    dt_plot[, total := shared + camdac_only]
    dt_plot[, prevalence :=
    belinksky_genes_prevalence[histology == hs][
      match(dt_plot$external_gene_name, gene), prevalence]
      ]
    dt_plot <- dt_plot[!is.na(prevalence)]
    dt_plot[, c("min", "max") := tstrsplit(prevalence, split = "-")]
    dt_plot[, c("min", "max") := lapply(.SD, function(x) as.numeric(x)/100),
          .SDcols = c("min", "max")]
    dt_plot[is.na(max), "max" := min]
    diffs <- dt_plot[,
      ifelse(abs(total - min) < abs(total - max),
             abs(total - min),
             abs(total - max))]
    diff_list[[length(diff_list) + 1]] <- diffs


  # save prevalence data separately
  dt_prevalence <- unique(dt_melt[, .(gene_name, hist, min, max)])

  # make plot
  p2 <- ggplot(dt_melt) +
      geom_col(aes(x = gene_name, y = value, color = variable, fill = variable),
      position = "stack", alpha = 0.25) +
      geom_point(data = dt_prevalence, aes(x = gene_name, y = min)) +
      geom_point(data = dt_prevalence, aes(x = gene_name, y = max)) +
      geom_linerange(data = dt_prevalence, aes(x = gene_name, ymin = min, ymax = max)) +
      scale_color_manual(name = "class",
        #values = c("bulk_only" = "lightgrey", "camdac_only" = "dodgerblue", "shared" = "black"),
        values = c("camdac_only" = "dodgerblue", "shared" = "black"),
        drop = FALSE) +
      scale_fill_manual(name = "class",
        values = c("camdac_only" = "dodgerblue", "shared" = "black"),
        drop = FALSE) +
      scale_y_continuous("patients (%)") +
      facet_wrap(. ~ hist) +
      theme_minimal() +
      coord_cartesian(ylim = c(0, 1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))

  if (hs != "NSCLC") p2 <- p2 + theme(legend.position = "none")
  assign(paste0("p2_", tolower(hs)), p2)

  # store data for summary metrics
  metrics_list[[length(metrics_list) + 1]] <- dt_plot
}

# save plot
gd <- arrangeGrob(p2_luad, p2_lusc, p2_nsclc, nrow = 1)
f_nm <- file.path("outputs", "DMR_recall_camdac_pp_bar_plot_tmp.")
ggsave(plot = gd, filename = paste0(f_nm, "svg"), device = "svg",
       width = 8, height = 2.75)
ggsave(plot = gd, filename = paste0(f_nm, "pdf"), device = "pdf",
       width = 8, height = 2.75)

# get the difference between this study and reported prevalence
all_diffs <- unlist(diff_list)
mean(all_diffs)
min(all_diffs)
max(all_diffs)

# get the mean, min and max camdac only fraction
all_data <- rbindlist(metrics_list)
print(all_data[, min(camdac_only)])
print(all_data[, max(camdac_only)])
print(all_data[, mean(camdac_only)])

# summarise across multiple samples by patient
dt_summary_by_gene[, num_regions := length(DMR),
  by = .(external_gene_name, method, patient, hist)]
dt_summary_by_patient <-
  dt_summary_by_gene[, .(
    m_diff_DMPs = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMR_ubiquity = ifelse(sum(!is.na(DMR)) > 0,
      sum(!is.na(DMR)) / num_regions[1], as.numeric(NA)),
    hyper_ubiquity =
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) > 0,
       sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) /
       num_regions[1], as.numeric(NA)),
        DMR = ifelse(sum(!is.na(DMR)) > 0, "DMR", NA_character_),
    DMR = ifelse(sum(!is.na(DMR)) > 0, "DMR", NA_character_),
    DMR_type =
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) > 0,
        "hyper",
      ifelse(sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE) > 0,
        "hypo",
      NA_character_))
      ),
  by = .(external_gene_name, method, patient, hist)]

# summarise DMRs by method across cluster of the same genes and patients
dt_summary_across_patient <-
  dt_summary_by_patient[, .(
    num_patients = length(patient),
    hyper_ubi_pp =
      sum(!is.na(hyper_ubiquity) & hyper_ubiquity == 1, na.rm = TRUE),
    DMRs_pp = sum(!is.na(DMR)),
    hyper_pp = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
    ),
by = .(external_gene_name, method, hist)]

# get the differences between bulk and purified data
histologies <- c("LUAD", "LUSC", "NSCLC")
metrics_list <- list()
for (hs in histologies) {
  # subset to histology of interest
  if (hs != "NSCLC") {
    dt_plot <- dt_summary_across_patient[hist == hs]
  } else {
    dt_plot <- dt_summary_across_patient
    dt_plot$hist <- "NSCLC"
  }

  # generate data from dot plot
  dt_plot <-
    dt_plot[,
      .(ubiq_all = sum(hyper_ubi_pp[method == "bulk_tumour"]) /
              sum(hyper_pp[method == "purified_tumour"]),
        ubiq_camdac_only = (sum(hyper_ubi_pp[method == "purified_tumour"])
              - sum(hyper_ubi_pp[method == "bulk_tumour"])) /
              sum(hyper_pp[method == "purified_tumour"])),
      by = .(external_gene_name, hist)]
  dt_plot[is.finite(ubiq_all) | is.finite(ubiq_camdac_only), "non_ubiq" :=
    (1 - ubiq_all - ubiq_camdac_only)]
  dt_plot[!is.finite(ubiq_all), ubiq_all := 0]
  dt_plot[!is.finite(ubiq_camdac_only), ubiq_camdac_only := 0]

  if (hs == "NSCLC") {
    print(dt_plot[, min(ubiq_camdac_only, na.rm = TRUE)])
    print(dt_plot[, max(ubiq_camdac_only, na.rm = TRUE)])
    print(dt_plot[, mean(ubiq_camdac_only, na.rm = TRUE)])
  }
  dt_melt <- melt(dt_plot, id.vars = c("external_gene_name", "hist"))

  # get expected prevalence
  dt_melt[, prevalence :=
    belinksky_genes_prevalence[histology == hs][
      match(dt_melt$external_gene_name, gene), prevalence]
      ]

  # remove genes without reported prevalence
  dt_melt <- dt_melt[!is.na(prevalence)]

  # order genes
  if (hs == "NSCLC") {
    gene_ord <- gene_ord_nsclc
  } else if (hs == "LUAD") {
    gene_ord <- gene_ord_luad
  } else {
    gene_ord <- gene_ord_lusc
  }
  dt_melt[, gene_name :=
    factor(external_gene_name, levels = gene_ord, ordered = TRUE)]

  # factorise colouring
  levs <- c("non_ubiq", "ubiq_camdac_only", "ubiq_all")
  dt_melt[, variable := factor(variable, levels = levs, ordered = TRUE)]

  # make plot
  p3 <- ggplot(dt_melt) +
      geom_col(aes(x = gene_name, y = value, color = variable, fill = variable),
      position = "stack", alpha = 0.5) +
      scale_color_manual(name = "class",
        values = c("non_ubiq" = "grey", "ubiq_camdac_only" = "#0CC40C", "ubiq_all" = "black"),
        drop = FALSE) +
      scale_fill_manual(name = "class",
        values = c("non_ubiq" = "grey", "ubiq_camdac_only" = "#0CC40C", "ubiq_all" = "black"),
        drop = FALSE) +
      scale_y_continuous("ubiquitous (%)") +
      facet_wrap(. ~ hist) +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))

  if (hs != "NSCLC") p3 <- p3 + theme(legend.position = "none")
  assign(paste0("p3_", tolower(hs)), p3)

  # store data for summary metrics
  metrics_list[[length(metrics_list) + 1]] <- dt_plot
}

# save plot
gd <- arrangeGrob(p3_luad, p3_lusc, p3_nsclc, nrow = 1, widths = c(1, 1, 1.15))
f_nm <- file.path("outputs", "DMR_recall_camdac_ubiquity_bar_plot.")
ggsave(plot = gd, filename = paste0(f_nm, "svg"), device = "svg",
       width = 8, height = 2.75)
ggsave(plot = gd, filename = paste0(f_nm, "pdf"), device = "pdf",
       width = 8, height = 2.75)

# get the mean, min and max camdac only fraction
all_data <- rbindlist(metrics_list)
print(all_data[, quantile(ubiq_camdac_only/(1 - non_ubiq), 0.25, na.rm=TRUE)])
print(all_data[, quantile(ubiq_camdac_only/(1 - non_ubiq), 0.75, na.rm=TRUE)])
print(all_data[, mean(ubiq_camdac_only/(1 - non_ubiq), na.rm=TRUE)])
print(all_data[, 1 - quantile(non_ubiq, 0.25, na.rm=TRUE)])
print(all_data[, 1 - quantile(non_ubiq, 0.75, na.rm=TRUE)])
print(all_data[, 1 - mean(non_ubiq, na.rm=TRUE)])

################################################################################
# new DMR discovery
################################################################################

# create results object
summary_list <- list()

# columns to select
keep_cols <-
  c("patient", "sample", "ID",
    "external_gene_name", "ensembl_gene_id",
    "method", "m_n", "m_n_DMPs", "m_diff_DMPs", "DMR",
    "CpG_counts", "DMP_counts", "consec_DMPs")

# iterate over bulk and purified tumuour DMR data
for (prefx in prefixs) {

  # set promoter DMR file name
  if (prefx == "purified_tumour") {
    f_nm <- file.path("outputs", "promoter_camdac_DMRs_nonas.rds")
    #f_nm <- file.path(tx_dir, "promoter_camdac_DMRs_nonas.rds")
  } else {
    f_nm <- file.path("outputs", "promoter_bulk_DMRs_nonas.rds")
    #f_nm <- file.path(tx_dir, "promoter_bulk_DMRs_nonas.rds")
  }

  # check if the processed file exists
  if (file.exists(f_nm)) {
    dt_prom_dmrs_final <- readRDS(f_nm)

    # compute the methylation difference
    if (prefx == "purified_tumour") {
      dt_prom_dmrs_final[, DMR := DMR_tn]
      dt_prom_dmrs_final[, m_diff_DMPs :=
        as.numeric(m_t_DMPs) - as.numeric(m_n_DMPs)]
    } else {
      dt_prom_dmrs_final[, DMR := DMR_bn]
      dt_prom_dmrs_final[, m_diff_DMPs :=
        as.numeric(m_b_DMPs) - as.numeric(m_n_DMPs)]
    }
    
    # annotate method
    dt_prom_dmrs_final$method <- prefx

    # add to list
    summary_list[[length(summary_list) + 1]] <-
      dt_prom_dmrs_final[, .SD, .SDcols = keep_cols]

  # otherwise, compile the gene-wise promoter methylation data
  } else {
    # load all genes
    load("inputs/hg19_external_gene_id_UCSC.RData")
    ucsc_gene_names <- genes$external_gene_name
    ucsc_gene_ids <- genes$ensembl_gene_id

    # set CAMDAC regions file name
    if (prefx == "purified_tumour") {
        f_nm2 <-
          file.path(tx_dir,
            "Groups_wise_CAMDAC_results_noNAs.fst")
    } else {
        f_nm2 <-
          file.path(tx_dir,
            "Groups_wise_bulk_tumour_CAMDAC_results_noNAs.fst")
    }
    dt_regions <- read.fst(path = f_nm2, as.data.table = TRUE)

    # select promoters DMRs
    dt_proms <- dt_regions[promoter == "promoter"]
    if(prefx == "purified_tumour"){
      dt_prom_dmrs <- dt_proms[!is.na(DMR_tn) & DMR_tn == "DMR"]
    } else {
      dt_prom_dmrs <- dt_proms[!is.na(DMR_bn) & DMR_bn == "DMR"]
    }

    # clean up
    rm(dt_proms, dt_regions)
    
    # get the promoter methylation
    lst <- mclapply(1:length(ucsc_gene_ids),
      function(i,
        genes = ucsc_gene_ids,
        dt_tmp = dt_prom_dmrs,
        names = ucsc_gene_names)
      for(gene in genes[i]) {
              if (which(genes == gene) %% 250 == 0) {
                cat(which(genes == gene) / 250, sep = "...")
              }
              tmp <- dt_tmp[grepl(gene, gene_prom_ids), ]
              if (nrow(tmp) > 0) {
                expr <- paste0("(^", gene, ",)|(;", gene, ",)")
                tmp <- tmp[grepl(expr, gene_prom_ids),]
              } else {
                tmp <- tmp[1,]
              }
              tmp$external_gene_name <- names[i]
              tmp$ensembl_gene_id <- genes[i]
              return(tmp)
    }, ucsc_gene_ids, dt_prom_dmrs, ucsc_gene_names, mc.cores = 10)
    for (k in length(lst):1) {
      if (is.na(lst[[k]]$patient[1])){
           lst[[k]] <- NULL
      }
    }

    # bind results
    dt_prom_dmrs_final <- rbindlist(lst)
    saveRDS(dt_prom_dmrs_final, file = f_nm)

    # compute the methylation difference
    if (prefx == "purified_tumour") {
      dt_prom_dmrs_final[, DMR := DMR_tn]
      dt_prom_dmrs_final[, m_diff_DMPs :=
        as.numeric(m_t_DMPs) - as.numeric(m_n_DMPs)]
    } else {
      dt_prom_dmrs_final[, DMR := DMR_bn]
      dt_prom_dmrs_final[, m_diff_DMPs :=
        as.numeric(m_b_DMPs) - as.numeric(m_n_DMPs)]
    }
    
    # annotate method
    dt_prom_dmrs_final$method <- prefx

    # add to list
    summary_list[[length(summary_list) + 1]] <-
      dt_prom_dmrs_final[, .SD, .SDcols = keep_cols]
  }
}

# bind the results
dt_proms <- rbindlist(summary_list)

# print the number of unique genes
length(unique(dt_proms$ensembl_gene_id))

# format histology, LTX and CRUKIDs
dt_proms[, hist := patient_info[match(dt_proms$patient, Patient), Histology]]
dt_proms[, patient_crukid :=
                patient_info[match(dt_proms$patient, Patient), CRUKID]]
dt_proms[, hist := ifelse(grepl("adeno", hist), "LUAD", "LUSC")]
dt_proms[, sample_ltxid := paste0(patient, "_", sample)]
dt_proms[, sample_crukid := paste0(patient_crukid, "_", sample)]

# note down DMR type
dt_proms[, DMR_type :=
  ifelse(sum(is.na(m_diff_DMPs)) == length(m_n), as.character(NA),
      ifelse(sum(m_diff_DMPs > 0, na.rm = TRUE) > 0,
      "hyper", "hypo")),
    by = .(patient, sample, ID, external_gene_name, ensembl_gene_id, method)]

# summarise DMRs by method by cluster across samples.
dt_summary_across_samples <-
  dt_proms[, .(
    mean_m_diff = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),

    # NSCLC per sample
    DMRs_ps = sum(!is.na(DMR)),
    hypo_ps = sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE),
    hyper_ps = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
    ),
  by = .(ID, external_gene_name, ensembl_gene_id, method, hist)]

# look at the genes most different between LUSC and LUAD
tmp <-
  dt_summary_across_samples[, .(
    total_dmrs = (DMRs_ps[hist == "LUAD"] + DMRs_ps[hist == "LUSC"]),
    total_hyper_luad_hypo_lusc =
      hyper_ps[hist == "LUAD"] + hypo_ps[hist == "LUSC"],
    total_hypo_luad_hyper_lusc =
      hypo_ps[hist == "LUAD"] + hyper_ps[hist == "LUSC"]
    ),
    by = .(ID, external_gene_name, ensembl_gene_id, method)]
tmp[total_dmrs > 100]
tmp[order(total_hyper_luad_hypo_lusc, decreasing = TRUE)][1:10]
tmp[order(total_hypo_luad_hyper_lusc, decreasing = TRUE)][1:10]

# summarise DMRs by method by cluster across patient
dt_summary_by_sample <-
  dt_proms[,
    .(DMR = ifelse(sum(!is.na(DMR)) == 0, NA_character_, "DMR"),
      DMR_type = ifelse(sum(!is.na(DMR)) == 0, NA_character_,
    ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) >
    sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE), "hyper",
    "hypo"))),
    by = .(patient, ID, external_gene_name, method, hist)]

# summary across patients
dt_summary_across_patients <-
  dt_summary_by_sample[, .(
    DMRs_pp = sum(!is.na(DMR)),
    hypo_pp = sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE),
    hyper_pp = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
  ),
  by = .(ID, external_gene_name, method, hist)]

# look at the genes most different between LUSC and LUAD
tmp <-
  dt_summary_across_patients[, .(
    total_dmrs = (DMRs_pp[hist == "LUAD"] + DMRs_pp[hist == "LUSC"]),
    total_hyper_all =
      hyper_pp[hist == "LUAD"] + hyper_pp[hist == "LUSC"],
    total_hypo_all =
      hypo_pp[hist == "LUAD"] + hypo_pp[hist == "LUSC"],
    total_hyper_luad_hypo_lusc =
      hyper_pp[hist == "LUAD"] + hypo_pp[hist == "LUSC"],
    total_hypo_luad_hyper_lusc =
      hypo_pp[hist == "LUAD"] + hyper_pp[hist == "LUSC"]
    ),
    by = .(ID, external_gene_name, method)]

tmp[total_dmrs > 35][method == "purified_tumour"]
tmp[total_hyper_all > 35][method == "purified_tumour"]
tmp[total_hyper_luad_hypo_lusc > 20 & total_hypo_all < 3 &
  total_hypo_luad_hyper_lusc < 3][method == "purified_tumour"][
    order(total_hyper_luad_hypo_lusc, decreasing = TRUE)]
dt_proms[external_gene_name == "KCTD1" & ID == "18_24130367_24131174",
  quantile(m_n_DMPs, na.rm = TRUE)]

# exatract top differing between hist
dt_hyper_luad_nodmr_lusc <-
  tmp[total_hyper_luad_hypo_lusc > 20 & total_hypo_all < 3 &
  total_hypo_luad_hyper_lusc < 3][
    order(total_hyper_luad_hypo_lusc, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
dt_nodmr_luad_hyper_lusc <-
  tmp[total_hyper_luad_hypo_lusc < 3 & total_hypo_all < 3 &
  total_hypo_luad_hyper_lusc > 12][
    order(total_hypo_luad_hyper_lusc, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
levs <- c(dt_nodmr_luad_hyper_lusc$external_gene_name,
        dt_hyper_luad_nodmr_lusc$external_gene_name)
gene_pattern <- paste0("(^", paste(levs, collapse = "$|^"), ")")
id_pattern <- paste0("(",
  paste(dt_nodmr_luad_hyper_lusc$ID, 
        dt_hyper_luad_nodmr_lusc$ID, collapse = "|", sep = "|"),
        ")")

# set plotting colors
cols <- c("hypo" = "blue2", "hyper" = "gold2") #, "mixed" = "olivedrab3")

# summarise across samples from the same patient
dt_plot <-
  dt_proms[grepl(gene_pattern, external_gene_name) & grepl(id_pattern, ID),
    .(DMR_type = ifelse(sum(!is.na(DMR)) == 0, NA_character_,
    ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) >
    sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE), "hyper",
    "hypo"))),
    by = .(patient_crukid, patient, ID, external_gene_name, method, hist)]
dt_plot[, gene_name := factor(external_gene_name, levels = levs, ordered = TRUE)]
dt_plot[, IDs := gsub("(CRUK00)|(CRUK000)", "", patient_crukid)]
dt_plot[, flag := DMR_type]
dt_plot[, flag2 := as.numeric(factor(flag, level = c("hyper", "hypo"), ordered = TRUE))]
f_nm <- file.path("outputs", paste0(prefx, "_DMR_by_hist_by_patient_normal_um.pdf"))
p1 <-
  ggplot(dt_plot[method == "purified_tumour"]) +
    geom_point(aes(x = IDs, y = gene_name, col = flag), shape = 15) +
    theme_classic() + theme(axis.text.x = element_text(angle=90, size=5)) +
    facet_grid(. ~ hist, scales = "free_x", space = "free") +
    scale_color_manual(name = "flag", values = cols, na.value = "white")
ggsave(plot = p1, filename = f_nm, device = "pdf", units = "in", width = 8, height = 2)
ggsave(plot = p1, filename = gsub(".pdf", ".svg", f_nm), device = "svg", units = "in", width = 8, height = 2)

# exatract top 10 differing between hist
dt_hyper_luad_hypo_lusc <-
  tmp[order(total_hyper_luad_hypo_lusc, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
dt_hypo_luad_hyper_lusc <-
  tmp[order(total_hypo_luad_hyper_lusc, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
levs <- c(dt_hypo_luad_hyper_lusc$external_gene_name,
        dt_hyper_luad_hypo_lusc$external_gene_name)
gene_pattern <- paste0("(^", paste(levs, collapse = "$|^"), ")")
id_pattern <- paste0("(",
  paste(dt_hypo_luad_hyper_lusc$ID, 
        dt_hyper_luad_hypo_lusc$ID, collapse = "|", sep = "|"),
        ")")

# summarise across samples from the same patient
dt_plot <-
  dt_proms[grepl(gene_pattern, external_gene_name) & grepl(id_pattern, ID),
    .(DMR_type = ifelse(sum(!is.na(DMR)) == 0, NA_character_,
    ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) >
    sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE), "hyper",
    "hypo"))),
    by = .(patient_crukid, patient, ID, external_gene_name, method, hist)]
dt_plot[, gene_name := factor(external_gene_name, levels = levs, ordered = TRUE)]
dt_plot[, IDs := gsub("(CRUK00)|(CRUK000)", "", patient_crukid)]
dt_plot[, flag := DMR_type]
dt_plot[, flag2 := as.numeric(factor(flag, level = c("hyper", "hypo"), ordered = TRUE))]
f_nm <- file.path("outputs", paste0(prefx, "_DMR_by_hist_by_patient_normal_im.pdf"))
p1 <-
  ggplot(dt_plot[method == "purified_tumour"]) +
    geom_point(aes(x = IDs, y = gene_name, col = flag), shape = 15) +
    theme_classic() + theme(axis.text.x = element_text(angle=90, size=5)) +
    facet_grid(. ~ hist, scales = "free_x", space = "free") +
    scale_color_manual(name = "flag", values = cols, na.value = "white")
ggsave(plot = p1, filename = f_nm, device = "pdf", units = "in", width = 8, height = 2)
ggsave(plot = p1, filename = gsub(".pdf", ".svg", f_nm), device = "svg", units = "in", width = 8, height = 2)

# exatract top 10 differing between hist
dt_hyper_all <-
  tmp[order(total_hyper_all, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
dt_hypo_all <-
  tmp[order(total_hypo_all, decreasing = TRUE)][
    method == "purified_tumour"][1:5,
    .(ID, external_gene_name)]
levs <- c(dt_hyper_all$external_gene_name, dt_hypo_all$external_gene_name)
gene_pattern <- paste0("(^", paste(levs, collapse = "$|^"), ")")
id_pattern <- paste0("(", paste(dt_hyper_all$ID, dt_hypo_all$ID,
        collapse = "|", sep = "|"),
        ")")

# summarise across samples from the same patient
dt_plot <-
  dt_proms[grepl(gene_pattern, external_gene_name) & grepl(id_pattern, ID),
    .(DMR_type = ifelse(sum(!is.na(DMR)) == 0, NA_character_,
    ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) >
    sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE), "hyper",
    "hypo"))),
    by = .(patient_crukid, patient, ID, external_gene_name, method, hist)]
dt_plot[, gene_name := factor(external_gene_name, levels = levs, ordered = TRUE)]
dt_plot[, IDs := gsub("(CRUK00)|(CRUK000)", "", patient_crukid)]
dt_plot[, flag := DMR_type]
dt_plot[, flag2 := as.numeric(factor(flag, level = c("hyper", "hypo"), ordered = TRUE))]
f_nm <- file.path("outputs", paste0(prefx, "_hyper_by_hist_by_patient.pdf"))
p1 <-
  ggplot(dt_plot[method == "purified_tumour"]) +
    geom_point(aes(x = IDs, y = gene_name, col = flag), shape = 15) +
    theme_classic() + theme(axis.text.x = element_text(angle=90, size=5)) +
    facet_grid(. ~ hist, scales = "free_x", space = "free") +
    scale_color_manual(name = "flag", values = cols, na.value = "white")
ggsave(plot = p1, filename = f_nm, device = "pdf", units = "in", width = 8, height = 2)
ggsave(plot = p1, filename = gsub(".pdf", ".svg", f_nm), device = "svg", units = "in",
  width = 8, height = 2)

# summarise across multiple regions by genes
dt_summary_by_gene <-
  dt_proms[, .(
    m_diff_DMPs = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMR = ifelse(sum(!is.na(DMR)) == 0, NA_character_, "DMR"),
    DMR_type = ifelse(sum(!is.na(DMR)) == 0, NA_character_,
      ifelse(sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE) >
      sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE), "hyper",
      "hypo")),
    CpG_counts = ifelse(sum(!is.na(DMR)) > 0,
      max(CpG_counts[!is.na(DMR)], na.rm = TRUE),
      max(CpG_counts, na.rm = TRUE)),
    DMP_counts = ifelse(sum(!is.na(DMR)) > 0,
      as.integer(max(DMP_counts[!is.na(DMR)], na.rm = TRUE)),
      NA_integer_),
    consec_DMPs = ifelse(sum(!is.na(DMR)) > 0,
      as.integer(max(consec_DMPs[!is.na(DMR)], na.rm = TRUE)),
      NA_integer_)
    ),
  by = .(external_gene_name, ensembl_gene_id, method, patient, sample, hist)]

# summarise DMRs by method across cluster of the same genes and samples.
dt_summary_across_samples <-
  dt_summary_by_gene[, .(
    mean_m_diff = mean(m_diff_DMPs[!is.na(DMR)], na.rm = TRUE),
    DMRs_ps = sum(!is.na(DMR)),
    hypo_ps = sum(!is.na(DMR) & DMR_type == "hypo", na.rm = TRUE),
    hyper_ps = sum(!is.na(DMR) & DMR_type == "hyper", na.rm = TRUE)
    ),
  by = .(external_gene_name, ensembl_gene_id, method, hist)]

################################################################################
# Per patient plot of key DMRs
################################################################################

# selected genes
genes=c("CDKN2A","PAX5","CHFR","MGMT","DAPK1",
        "RASSF1","CDH1","CDH13","TIMP3",
        "LAMA3","LAMB3","LAMC2")
# NORE1 is excluded due to not being covered by RRBS 
# CASP8 0% NSCLC, "FAS" and "TRAIL1" SCLC
# ref: Belinsky, Nature Reviews Cancer, volume 4, 707–717 (2004)
# + GATA4, GATA5, PAX5α and PAX5β # 10.1158/1940-6207.CAPR-17-0177
# + APC, ESR1, GSTP1, HIC1, IGFBP3, NORE1A, PYCARD, RARB2, # 10.1016/B978-0-12-380866-0.60002-2
# TP73

for (prefx in prefixs) {
    
    # load DMRs
    if (prefx == "purified_tumour") {
        f_nm2 <- paste0(tx_dir, "Groups_wise_CAMDAC_results_noNAs.fst")
    } else {
        f_nm2 <- paste0(tx_dir, "Groups_wise_bulk_tumour_CAMDAC_results_noNAs.fst")  
    }
    df_final <- read.fst(path = f_nm2, as.data.table = TRUE)

    # select promoters DMRs
    df_proms <- df_final[promoter == "promoter",] 
    if (prefx == "purified_tumour") {
        df_proms <-  df_proms[!is.na(DMR_tn) & DMR_tn == "DMR",]  
        #df_proms <- df_proms[!is.na(DMR)&DMR=="DMR",] 
    } else {
        df_proms <-  df_proms[!is.na(DMR_bn) & DMR_bn == "DMR",]  
    }

    # subset patients 
    df_proms[, PXID := paste(patient, sample , sep="_")]

    # assign hist
    df_proms[, hist := patient_info$Histology[match(PXID, patient_info$LTX_ID)]]
    df_proms[, hist := ifelse(grepl("quamous", hist), "LUSC", "LUAD")]
  
    results_list <- list()
    genes_list <- list()
    for (gene in genes) {
        df <- data.table(name = gene, num_samples = 0, num_patient = 0)
        expr <- paste0("(^",gene,",)|(;",gene,",)")
        tmp <- df_proms[grepl(expr, gene_prom_ids), ]
  
        if (prefx == "purified_tumour" & nrow(tmp) > 0){
        df_gene <- data.table(name = gene, 
                            chr = unique(tmp$chromosome), 
                            min = min(tmp$start), 
                            max = max(tmp$start))
        genes_list[[length(genes_list)+1]] <- df_gene
    }
    
    # check that DMR is hyper and not hypo
    # if m_diff > 0 in m_t-m_n is is also >0 in m_b - m_n
    if (prefx == "purified_tumour") {
      tmp <- tmp[grepl("hyper", DMP_tn), ] 
    } else {
      tmp <- tmp[grepl("hyper", DMR_type), ]
    }
    
    df$num_samples <- tmp[, length(unique(PXID))]
    df$num_samples_LUAD <- tmp[hist=="LUAD", length(unique(PXID))]
    df$num_samples_LUSC <- tmp[hist=="LUSC", length(unique(PXID))]
    
    df$num_patient <- tmp[, length(unique(patient))]
    df$num_patient_LUAD <- tmp[hist=="LUAD", length(unique(patient))]
    df$num_patient_LUSC <- tmp[hist=="LUSC", length(unique(patient))]
    
    results_list[[length(results_list) + 1]] <- df
  }
  df1 <- rbindlist(results_list)
  assign(paste0("df1_", prefx), df1)
  print(df1)
  
  if (prefx == "purified_tumour") {
    df_genes <- rbindlist(genes_list)
    rm(df_gene, genes_list)
  }
  
  results_list <- list()
  for(gene in genes){
    df <- data.frame(name = gene)
    df[, 2:39] <- 0
    colnames(df)[2:39] <- unique(patient_info$Patient)
    
    expr <- paste0("(^",gene,",)|(;",gene,",)")
    tmp <- df_proms[grepl(expr, gene_prom_ids), ]
    if(prefx == "purified_tumour"){
      tmp <- tmp[grepl("hyper", DMP_tn), ] 
    } else {
      tmp <- tmp[grepl("hyper", DMR_type), ]
    }
    tmp <- tmp[, .(l = length(unique(PXID))), .(patient)]
    
    tmp2 <- patient_info[, .(l2 = length(unique(LTX_ID))), .(Patient)]
    tmp2$l <- tmp$l[match(tmp2$Patient,tmp$patient)]
    tmp2[, ratio := ifelse(is.na(l2), 0, l/l2)]
    
    df[, 2:39] <- tmp2$ratio
    results_list[[length(results_list)+1]] <- df
  }
  df2 <- rbindlist(results_list)
  assign(paste0("df2_", prefx), df2)
  mean(colSums(!is.na(df2[, 2:39]), na.rm=TRUE))
  table(colSums(!is.na(df2[, 2:39]), na.rm=TRUE))

  rm(df, df1, df2, gene, tmp, results_list, tmp2)
}


for(patient_id in patient_ids){
    path=paste0("~/rdeocd/txMeth/",patient_id,"/")
    path_patient=path
    
    files <- list.dirs(paste0("~/rdeocd/txMeth/",patient_id,"/archive/Allelecounts"),recursive = FALSE)
    sample_ids <- unlist(lapply(1:length(files), function(i, files) strsplit(files, split="[//]")[[i]][9], files));rm(files)
    sample_ids <- sample_ids[grepl(pattern="^R[0-9]", sample_ids)]

      for(sample_id in sample_ids){
        print(patient_id)
        print(sample_id)
        
        # set output path
        f_nm <- paste0(path_patient, "/Methylation/",sample_id, "/re_segmented/purified_tumour_DMP_calls.RData")

        # check files exists
        if(file.exists(f_nm)){
          path_output <- paste0(path_patient, "Methylation/",sample_id,"/re_segmented/")
          load(f_nm)
        }
        if(!file.exists(f_nm)){
          f_nm <- paste0(path_patient,"/Methylation/",sample_id,"/purified_tumour_DMP_calls.RData")
          path_output <-paste0(path_patient, "Methylation/",sample_id,"/")
          load(f_nm)
        }
        
        unlink(paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/prom_plot_DMRs/*.pdf"), recursive = T)
        suppressWarnings(dir.create(paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/prom_plot_DMRs/"), recursive = T))
        genes <- df_genes$name[is.finite(df_genes$min)]
        df_mt <- data.table(df_mt)
        df_mt[, m_t := m_t_corr] 
        
        for(gene in genes){
          tmp <- df_genes[name==gene,]
          # df_mt_sub <- df_proms[end>=tmp$min & start<=tmp$max &chromosome==tmp$chr & 
          #                         patient==patient_id & sample==sample_id,]
          df_mt_sub <- df_mt[end>=tmp$min & start<=tmp$max &chromosome==tmp$chr,]
          df_mt_sub$DMP_b <- ifelse(is.na(df_mt_sub$DMP_b), "NA", df_mt_sub$DMP_b)
          df_mt_sub$DMP_t <- ifelse(is.na(df_mt_sub$DMP_t), "NA", df_mt_sub$DMP_t)
        
          df_mt_sub$dist <- c(0,diff(df_mt_sub$start))
          flag <- ifelse((df_mt_sub$dist > 250), FALSE, TRUE)
          flag[1] <- F
          df_mt_sub$cluster_id <- as.character(cumsum(flag==FALSE))
          rm(flag)
        
        for(id in unique(df_mt_sub$cluster_id)){
          df_mt_sub_sub <- df_mt_sub[cluster_id==id, ]
        if(nrow(df_mt_sub_sub)>4){
          fc = factor(c("N","B","T","hyper","hypo","NA"),levels=c("N","B","T","hyper","hypo","NA"),ordered=TRUE)
          ggplot(df_mt_sub_sub)+geom_point(aes(x=start,y=m_n,col="N"),size=3)+
              geom_point(aes(x=start,y=m_t,colour=DMP_t),size=3)+
              geom_point(aes(x=start,y=m_b,colour=DMP_b),size=3)+theme_classic()+
              scale_colour_manual(name = "Methylation\nRate", 
                                  values = c("N"="rosybrown","B"="mediumpurple3","T"="lightgrey","hyper"="chartreuse3","hypo"="cornflowerblue",
                                             "NA"="lightgrey"), labels = c("N"="Normal","B"="Bulk","T"="Tumour","hyper"="Tumour\nhyper-\nmethylated",
                                                                           "hypo"="Tumour\nhypo-\nmethylated", "NA"="Tumour non-\ndifferentially\nmethylated"),breaks = fc)+
              # geom_text(data=tmp,aes(x=min+((min+end)/2),y=1.1,label=paste0(gene, " promoter")),col="black",size=3,hjust=0)+
              # geom_label(data=tmp,aes(x=min+0.2*((min+end)/2),y=1.1,label=paste0(patient_id, "-", sample_id)),col="black", size=3)+
              geom_ribbon(aes(x=start, ymin=ifelse(m_t_low<0,0,m_t_low),ymax=ifelse(m_t_high>1,1,m_t_high)),alpha=0.5,fill="lightgrey")+
              geom_ribbon(aes(x=start, ymin=ifelse(m_n_low<0,0,m_n_low),ymax=ifelse(m_n_high>1,1,m_n_high)),alpha=0.5,fill="rosybrown")+
              geom_ribbon(aes(x=start, ymin=ifelse(m_b_low<0,0,m_b_low),ymax=ifelse(m_b_high>1,1,m_b_high)),alpha=0.5, fill="mediumpurple3")+
              geom_line(aes(x=start,y=m_t,col="T"),alpha=0.75)+xlab("Chromosome Position (bp)")+
              geom_line(aes(x=start,y=m_b, col="B"),alpha=0.75)+geom_line(aes(x=start,y=m_n, col="N"))+
              scale_y_continuous(name = "methylation rate", breaks = seq(0,1,0.1), limits = c(-0.1,1.1)) +
              ggsave(filename = paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/prom_plot_DMRs/", 
                                       gene,"_prom_plot_",id,".pdf"),device = "pdf",width = 8, height = 4)
          
          fc = factor(c("B","T","hyper","hypo","NA"),levels=c("B","T","hyper","hypo","NA"),ordered=TRUE)
          ggplot(df_mt_sub_sub)+
            geom_point(aes(x=start,y=m_t-m_n,colour=DMP_t),size=3)+
            geom_point(aes(x=start,y=m_b-m_n,colour=DMP_b),size=3)+theme_classic()+
            scale_colour_manual(name = "Methylation\nRate", 
                                values = c("B"="mediumpurple3","T"="lightgrey","hyper"="chartreuse3","hypo"="cornflowerblue", "NA"="lightgrey"), 
                                labels = c("B"="Bulk","T"="Tumour","hyper"="Tumour\nhyper-\nmethylated","hypo"="Tumour\nhypo-\nmethylated", 
                                           "NA"="Tumour non-\ndifferentially\nmethylated"), breaks = fc)+
            # geom_text(data=tmp,aes(x=min+((min+end)/2),y=1.1,label=paste0(gene, " promoter")),col="black",size=3,hjust=0)+
            # geom_label(data=tmp,aes(x=min+0.2*((min+end)/2),y=1.1,label=paste0(patient_id, "-", sample_id)),col="black", size=3)+
            geom_ribbon(aes(x=start, ymin=ifelse(m_t_low-m_n<(-1),-1,m_t_low-m_n),ymax=ifelse(m_t_high-m_n>1,1,m_t_high-m_n)),alpha=0.5,fill="lightgrey")+
            geom_ribbon(aes(x=start, ymin=ifelse(m_b_low-m_n<(-1),-1,m_b_low-m_n),ymax=ifelse(m_b_high-m_n>1,1,m_b_high-m_n)),alpha=0.5, fill="mediumpurple3")+
            geom_ribbon(aes(x=start, ymin=(0-(m_n-m_n_low)),ymax=(m_n_high-m_n)),alpha=0.5, fill="rosybrown")+
            geom_line(aes(x=start,y=m_t-m_n,col="T"),alpha=0.75)+geom_line(aes(x=start,y=m_b-m_n,col="T"),alpha=0.75)+
            xlab("Chromosome Position (bp)")+
            scale_y_continuous(name = "methylation rate", breaks = seq(-1,1,0.2), limits = c(-1,1.1)) +
            ggsave(filename = paste0("~/rdeocd/txMeth/",patient_id,"/Methylation/",sample_id,"/prom_plot_DMRs/m_diff_", 
                                     gene,"_prom_plot_",id,".pdf"),device = "pdf",width = 8, height = 4)
        }
      }
    }
  }
}
rm(list=ls(pattern="(df)|(patient)"))
rm(tmp,fc,gene,genes)

# END