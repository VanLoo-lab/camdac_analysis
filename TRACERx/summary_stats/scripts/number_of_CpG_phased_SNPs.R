#!/usr/bin/env Rscript

#' @description Investigating the number of CpG sites with phased SNPs
#' in the CAMDAC dataset
#' This script loads all data and counts:
#' 1. Unique SNPs (CHR, POS_genotype)
#' 2. Unique CpGs (CHR, start, end) where total_counts_m_n_2 > 0 and total_counts_m_n_1 > 0
#'
#' @author Lili Larose Cadieux

# Set up error handling and logging
options(error = function() {
  cat("Error occurred. Check the log file for details.\n")
  sink()
  q(status = 1)
})

# Redirect stdout to log file
if (!dir.exists("outputs")) {
  dir.create("outputs", recursive = TRUE)
}
sink("outputs/log.out", append = FALSE, split = TRUE)

# Load library
x <- c("data.table", "stringr", "ggplot2", "gridExtra", "dplyr")

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

# Set base directory
base_dir <- "/home/lili/rdeocd/txMeth"

# Function to get all patient directories
get_patient_dirs <- function(base_dir) {
  patient_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  patient_dirs <- patient_dirs[grepl("LTX", basename(patient_dirs))]
  return(patient_dirs)
}

# Function to get all sample directories for a patient
get_sample_dirs <- function(patient_dir) {
  asm_dir <- file.path(patient_dir, "Allele_specific_analysis")
  if (!dir.exists(asm_dir)) return(NULL)
  
  sample_dirs <- list.dirs(asm_dir, full.names = TRUE, recursive = FALSE)
  sample_dirs <- sample_dirs[basename(sample_dirs) %in% c("N", paste0("R", 1:8))]
  return(sample_dirs)
}

# Function to load and analyze allele-specific methylation data for a sample
load_and_analyze_asm_data <- function(sample_dir) {
  patient <- basename(dirname(dirname(sample_dir)))
  sample <- basename(sample_dir)
  
  # Look for the RData file
  rdata_file <- file.path(sample_dir,
    paste0("allele_specific_methylation_rates_", sample, ".RData"))
  
  if (!file.exists(rdata_file)) {
    return(NULL)
  }
  
  tryCatch({
    # Load the data
    load(rdata_file)
    
    # The object name is consistently 'df_phased_methylation' across all files
    if (exists("df_phased_methylation")) {
      data <- df_phased_methylation
      
      # Analyze unique counts directly
      counts <- analyze_unique_counts(data)
      
      # Create result list with counts and metadata
      result <- list(
        patient = patient,
        sample = sample,
        sample_id = paste0(patient, "_", sample),
        unique_snps = counts$snps,
        unique_cpgs = counts$cpgs,
        total_rows = nrow(data)
      )
      
      # Clear the data from memory
      rm(data, df_phased_methylation)
      gc()
      
      return(result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    return(NULL)
  })
}

# Function to analyze unique SNPs and CpGs
analyze_unique_counts <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    return(list(snps = 0, cpgs = 0))
  }
  
  # Count unique SNPs (CHR, POS_genotype)
  unique_snps <- data %>%
    select(CHR, POS_genotype) %>%
    distinct() %>%
    nrow()
  
  # Count unique CpGs (CHR, start, end) where total_counts_m_n_2 > 0 and total_counts_m_n_1 > 0
  unique_cpgs <- data %>%
    filter(total_counts_m_n_2 > 0 & total_counts_m_n_1 > 0 & 
      !is.na(total_counts_m_n_2) & !is.na(total_counts_m_n_1)) %>%
    select(CHR, start, end) %>%
    distinct() %>%
    nrow()
  
  return(list(snps = unique_snps, cpgs = unique_cpgs))
}

# Main execution
cat("Loading and analyzing allele-specific methylation data...\n")

# Get all patient directories
patient_dirs <- get_patient_dirs(base_dir)
cat("Found", length(patient_dirs), "patient directories\n")

# Initialize results storage as named list
all_results <- list()
total_unique_snps <- 0
total_unique_cpgs <- 0
patient_count <- 0
sample_count <- 0

# Process each patient
for (patient_dir in patient_dirs) {
  patient <- basename(patient_dir)
  cat("\nProcessing patient:", patient, "\n")
  
  # Get sample directories
  sample_dirs <- get_sample_dirs(patient_dir)
  
  if (is.null(sample_dirs) || length(sample_dirs) == 0) {
    cat("  No sample directories found\n")
    next
  }
  
  cat("  Found", length(sample_dirs), "sample directories\n")
  
  # Process each sample
  for (sample_dir in sample_dirs) {
    sample <- basename(sample_dir)
    cat("  Loading sample:", sample, "... ")
    
    # Load and analyze the data
    result <- load_and_analyze_asm_data(sample_dir)
    
    if (!is.null(result)) {
      # Store result in named list
      all_results[[result$sample_id]] <- result
      
      total_unique_snps <- total_unique_snps + result$unique_snps
      total_unique_cpgs <- total_unique_cpgs + result$unique_cpgs
      sample_count <- sample_count + 1
      cat("SUCCESS (", result$unique_snps, "SNPs,", result$unique_cpgs, "CpGs)\n")
    } else {
      cat("FAILED\n")
    }
  }
  
  patient_count <- patient_count + 1
}

# Summary
cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("SUMMARY:\n")
cat("Patients processed:", patient_count, "\n")
cat("Samples loaded:", sample_count, "\n")
cat("Total unique SNPs across all samples:", total_unique_snps, "\n")
cat("Total unique CpGs across all samples:", total_unique_cpgs, "\n")

# Show detailed results
if (length(all_results) > 0) {
  cat("\nDetailed results by sample:\n")
  cat("Sample_ID\t\tUnique_SNPs\tUnique_CpGs\tTotal_Rows\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  for (sample_id in names(all_results)) {
    result <- all_results[[sample_id]]
    cat(sprintf("%-20s\t%d\t\t%d\t\t%d\n", 
                sample_id, 
                result$unique_snps, 
                result$unique_cpgs, 
                result$total_rows))
  }
}

# Save results
cat("\nSaving results...\n")
save(all_results, file = "outputs/unique_snps_cpgs_analysis.RData")

# Create summary data frame
results_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    patient = x$patient,
    sample = x$sample,
    unique_snps = x$unique_snps,
    unique_cpgs = x$unique_cpgs,
    total_rows = x$total_rows,
    stringsAsFactors = FALSE
  )
}))

write.csv(results_df, file = "outputs/unique_snps_cpgs_summary.csv", row.names = FALSE)

# Also save as TSV file
write.table(results_df, file = "outputs/unique_snps_cpgs_summary.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Results saved to:\n")
cat("- unique_snps_cpgs_analysis.RData (R object)\n")
cat("- unique_snps_cpgs_summary.csv (CSV file)\n")
cat("- unique_snps_cpgs_summary.tsv (TSV file)\n")

# Compute the number of SNPs per patient and CpG per patient
results_df %>%
  filter(sample == "N") %>%
  summarise(total_snps = mean(unique_snps),
            total_cpgs = mean(unique_cpgs)) %>%
  print()

cat("\nAnalysis completed successfully!\n") 

# Close the log file and exit cleanly
sink()
cat("Analysis completed. Check outputs/log.out for detailed results.\n")
q(status = 0)