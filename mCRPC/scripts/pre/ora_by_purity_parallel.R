#!/camp/apps/eb/software/R/3.6.2-foss-2019b/bin/Rscript
#SBATCH --job-name=recc_dmrs
#SBATCH --output=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/logs/ora/241006_%A_%a.out
#SBATCH --error=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/logs/ora/241006_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Aim: Over-enrichment analysis on bulk and pure promoter DMRs for mCRPC, LUAD, and LUSC samples
# ORA analysis tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html
#  This script generates the ORA by purity output, but takes a parameter to run one of the sample lists through its pipeline

# Get the slurm array task id as a variable
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Setup -----
source("src/fora.R")
library_ora()
library(logger)
library(ggplot2)

# Files and directories
log_info("Running ORA by purity")
outdir <- fs::dir_create("data/pre/gsea_500")
wgbs_root <- "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/dmrs/"
rrbs_root <- "/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMRs/"
wgbs_pp <- fread("/camp/home/mensahn/projects/paper_camdac/data/out/qc/meta_mcrpc.tsv") %>% filter(tPASS == T) # Purity and ploidy file
rrbs_pp <- fread("~/projects/paper_camdac/data/out/qc/rrbs_pp.tsv") %>% dplyr::mutate(id = paste0(tracerx_id, "_", Region_ID))
luad_pp <- rrbs_pp %>% dplyr::filter(id %in% get_rrbs_subtype('LUAD'))
lusc_pp <- rrbs_pp %>% dplyr::filter(id %in% get_rrbs_subtype('LUSC'))

log_info("Generating lists of samples for each threshold")

# Threshold purity by the median for each group only (0.5 quantile)
w_lq <- wgbs_pp %>%
  filter(purity <= quantile(purity, 0.5)) %>%
  pull(sample_id)
w_hq <- wgbs_pp %>%
  filter(purity > quantile(purity, 0.5)) %>%
  pull(sample_id)
rluad_lq <- luad_pp %>%
  filter(rrbs_rho <= quantile(rrbs_rho, 0.5)) %>%
  pull(id)
rluad_hq <- luad_pp %>%
  filter(rrbs_rho> quantile(rrbs_rho, 0.5)) %>%
  pull(id)
rlusc_lq <- lusc_pp %>%
  filter(rrbs_rho <= quantile(rrbs_rho, 0.5)) %>%
  pull(id)
rlusc_hq <- lusc_pp %>%
  filter(rrbs_rho> quantile(rrbs_rho, 0.5)) %>%
  pull(id)


# Load all WGBS and RRBS (luad/lusc) files
log_info("Loading all WGBS and RRBS files")
wgbs_root <- "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/dmrs/"
rrbs_root <- "/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMRs/"
wgbs_pure <- fs::dir_ls(wgbs_root, recurse = T, glob = "*pure_dmr.fst")
wgbs_bulk <- fs::dir_ls(wgbs_root, recurse = T, glob = "*bulk_dmr.fst")
lusc_pure <- fs::dir_ls(rrbs_root, recurse = T, glob = "*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt = "LUSC")
lusc_bulk <- fs::dir_ls(rrbs_root, recurse = T, glob = "*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt = "LUSC")
luad_pure <- fs::dir_ls(rrbs_root, recurse = T, glob = "*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt = "LUAD")
luad_bulk <- fs::dir_ls(rrbs_root, recurse = T, glob = "*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt = "LUAD")


# Lists: Generate lists for WGBS, luad and lusk by filtering r_* file lists

## WGBS file lists
wgbs_pure_lq <- filter_samples(wgbs_pure, w_lq)
wgbs_pure_hq <- filter_samples(wgbs_pure, w_hq)
wgbs_bulk_lq <- filter_samples(wgbs_bulk, w_lq)
wgbs_bulk_hq <- filter_samples(wgbs_bulk, w_hq)

## LUAD file lists
luad_pure_lq <- filter_tx_samples(luad_pure, rluad_lq)
luad_pure_hq <- filter_tx_samples(luad_pure, rluad_hq)
luad_bulk_lq <- filter_tx_samples(luad_bulk, rluad_lq)
luad_bulk_hq <- filter_tx_samples(luad_bulk, rluad_hq)

## LUSC file lists
lusc_pure_lq <- filter_tx_samples(lusc_pure, rlusc_lq)
lusc_pure_hq <- filter_tx_samples(lusc_pure, rlusc_hq)
lusc_bulk_lq <- filter_tx_samples(lusc_bulk, rlusc_lq)
lusc_bulk_hq <- filter_tx_samples(lusc_bulk, rlusc_hq)

# Count the number of samples in each category and save based on pure files only
ccount <- data.frame(
  category = c( "w_lq", "w_hq",  "luad_lq", "luad_hq", "lusc_lq", "lusc_hq"),
  count = c(length(wgbs_pure_lq), length(wgbs_pure_hq),length(luad_pure_lq), length(luad_pure_hq),  length(lusc_pure_lq), length(lusc_pure_hq))
)
data.table::fwrite(ccount, fs::path(outdir, "sample_counts.tsv"), sep = "\t", quote = F)

log_info("Filtering files to sample lists")

## Make cache file listsings
wpure_hq_out <- fs::path(outdir, "wgbs_pure_hq.qs")
wpure_lq_out <- fs::path(outdir, "wgbs_pure_lq.qs")
wbulk_hq_out <- fs::path(outdir, "wgbs_bulk_hq.qs")
wbulk_lq_out <- fs::path(outdir, "wgbs_bulk_lq.qs")
luad_pure_hq_out <- fs::path(outdir, "luad_pure_hq.qs")
luad_pure_lq_out <- fs::path(outdir, "luad_pure_lq.qs")
luad_bulk_hq_out <- fs::path(outdir, "luad_bulk_hq.qs")
luad_bulk_lq_out <- fs::path(outdir, "luad_bulk_lq.qs")
lusc_pure_hq_out <- fs::path(outdir, "lusc_pure_hq.qs")
lusc_pure_lq_out <- fs::path(outdir, "lusc_pure_lq.qs")
lusc_bulk_lq_out <- fs::path(outdir, "lusc_bulk_lq.qs")
lusc_bulk_hq_out <- fs::path(outdir, "lusc_bulk_hq.qs")

## Generate a selector for the analysis parameters based on the task ID. It should give a vector that can be passed to the pipe
task_selector <- list(
  "1" = list(wgbs_pure_lq, wgbs_bulk_lq, wpure_lq_out, wbulk_lq_out, fs::path(outdir, "ORA_wgbs_lq.qs")),
  "2" = list(wgbs_pure_hq, wgbs_bulk_hq, wpure_hq_out, wbulk_hq_out, fs::path(outdir, "ORA_wgbs_hq.qs")),
  "3" = list(luad_pure_lq, luad_bulk_lq, luad_pure_lq_out, luad_bulk_lq_out, fs::path(outdir, "ORA_luad_lq.qs")),
  "4" = list(luad_pure_hq, luad_bulk_hq, luad_pure_hq_out, luad_bulk_hq_out, fs::path(outdir, "ORA_luad_hq.qs")),
  "5" = list(lusc_pure_lq, lusc_bulk_lq, lusc_pure_lq_out, lusc_bulk_lq_out, fs::path(outdir, "ORA_lusc_lq.qs")),
  "6" = list(lusc_pure_hq, lusc_bulk_hq, lusc_pure_hq_out, lusc_bulk_hq_out, fs::path(outdir, "ORA_lusc_hq.qs"))
)
TASK_PARAMS <- task_selector[[task_id]]

# Log and quit if the ORA output file exists
if (fs::file_exists(TASK_PARAMS[[5]])) {
  log_info("ORA output file exists for {TASK_PARAMS[[5]]}. Exiting")
  quit()
}

# Select the function, one of pipe_ora_wgbs or pipe_ora_rrbs based on whether
function_selector <- ifelse(stringr::str_detect(TASK_PARAMS[[5]], "wgbs"), pipe_ora_wgbs, pipe_ora_rrbs)

# Run ORA for the selected task
log_info("Running ORA for {task_id}")
log_info("Params = {str(TASK_PARAMS,1)}")
ora_res <- function_selector(TASK_PARAMS[[1]], TASK_PARAMS[[2]], TASK_PARAMS[[3]], TASK_PARAMS[[4]], TASK_PARAMS[[5]])
log_info("ORA complete")
