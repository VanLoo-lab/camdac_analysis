
# Aim: Over-enrichment analysis on bulk and pure promoter DMRs for mCRPC, LUAD, and LUSC samples
# ORA analysis tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html

# Setup -----
source('src/fora.R')
library_ora()
library(logger)
library(ggplot2)

# Files and directories
log_info("Running ORA by purity")
outdir = fs::dir_create("data/pre/gsea")
wgbs_root = "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/dmrs/"
rrbs_root = "/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMRs/"
wgbs_pp = fread('/camp/home/mensahn/projects/paper_camdac/data/out/qc/meta_mcrpc.tsv') %>% filter(tPASS==T)# Purity and ploidy file
rrbs_pp = fread('~/projects/paper_camdac/data/out/qc/rrbs_pp.tsv') %>% dplyr::mutate(id=paste0(tracerx_id, "_", Region_ID))

log_info("Generating lists of samples for each threshold")
# WGBS: Generate lists of bulk and pure sample filters for each purity threshold:
# hf:>0.5 and lf:<=0.5; and low:0-.33, mid:>.33-.66, high:>.66
w_lf = wgbs_pp %>% filter(purity <= .5) %>% pull(sample_id)
w_hf = wgbs_pp %>% filter(purity > .5) %>% pull(sample_id)
w_low = wgbs_pp %>% filter(purity <= .33) %>% pull(sample_id)
w_mid = wgbs_pp %>% filter(purity > .33 & purity <= .66) %>% pull(sample_id)
w_high = wgbs_pp %>% filter(purity > .66) %>% pull(sample_id)

# RRBS: Generate lists of bulk and pure sample filters for each purity threshold
# hf:>0.5 and lf:<=0.5; and low:0-.33, mid:>.33-.66, high:>.66
r_lf = rrbs_pp %>% filter(rrbs_rho <= .5) %>% pull(id)
r_hf = rrbs_pp %>% filter(rrbs_rho > .5) %>% pull(id)
r_low = rrbs_pp %>% filter(rrbs_rho <= .33) %>% pull(id)
r_mid = rrbs_pp %>% filter(rrbs_rho > .33 & rrbs_rho <= .66) %>% pull(id)
r_high = rrbs_pp %>% filter(rrbs_rho > .66) %>% pull(id)

# Count the number of samples in each category and save in atable
ccount = data.frame(
    category=c('w_lf', 'w_hf', 'w_low', 'w_mid', 'w_high', 'r_lf', 'r_hf', 'r_low', 'r_mid', 'r_high'),
    count=c(length(w_lf), length(w_hf), length(w_low), length(w_mid), length(w_high), length(r_lf), length(r_hf), length(r_low), length(r_mid), length(r_high))
)
data.table::fwrite(ccount, fs::path(outdir, 'sample_counts.tsv'), sep='\t', quote=F)

# Load all WGBS and RRBS (luad/lusc) files
log_info("Loading all WGBS and RRBS files")
wgbs_root = "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/dmrs/"
rrbs_root = "/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMRs/"
wgbs_pure = fs::dir_ls(wgbs_root, recurse=T, glob="*pure_dmr.fst")
wgbs_bulk= fs::dir_ls(wgbs_root, recurse=T, glob="*bulk_dmr.fst")
lusc_pure = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUSC")
lusc_bulk = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUSC")
luad_pure = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUAD")
luad_bulk = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUAD")

log_info("Filtering files to sample lists")
## WGBS file lists
wgbs_pure_lf = filter_samples(wgbs_pure, w_lf)
wgbs_pure_hf = filter_samples(wgbs_pure, w_hf)
wgbs_pure_low = filter_samples(wgbs_pure, w_low)
wgbs_pure_mid = filter_samples(wgbs_pure, w_mid)
wgbs_pure_high = filter_samples(wgbs_pure, w_high)
wgbs_bulk_lf = filter_samples(wgbs_bulk, w_lf)
wgbs_bulk_hf = filter_samples(wgbs_bulk, w_hf)
wgbs_bulk_low = filter_samples(wgbs_bulk, w_low)
wgbs_bulk_mid = filter_samples(wgbs_bulk, w_mid)
wgbs_bulk_high = filter_samples(wgbs_bulk, w_high)

## LUAD file lists
luad_pure_lf = filter_tx_samples(luad_pure, r_lf)
luad_pure_hf = filter_tx_samples(luad_pure, r_hf)
luad_pure_low = filter_tx_samples(luad_pure, r_low)
luad_pure_mid = filter_tx_samples(luad_pure, r_mid)
luad_pure_high = filter_tx_samples(luad_pure, r_high)
luad_bulk_lf = filter_tx_samples(luad_bulk, r_lf)
luad_bulk_hf = filter_tx_samples(luad_bulk, r_hf)
luad_bulk_low = filter_tx_samples(luad_bulk, r_low)
luad_bulk_mid = filter_tx_samples(luad_bulk, r_mid)
luad_bulk_high = filter_tx_samples(luad_bulk, r_high)

## LUSC file lists
lusc_pure_lf = filter_tx_samples(lusc_pure, r_lf)
lusc_pure_hf = filter_tx_samples(lusc_pure, r_hf)
lusc_pure_low = filter_tx_samples(lusc_pure, r_low)
lusc_pure_mid = filter_tx_samples(lusc_pure, r_mid)
lusc_pure_high = filter_tx_samples(lusc_pure, r_high)
lusc_bulk_lf = filter_tx_samples(lusc_bulk, r_lf)
lusc_bulk_hf = filter_tx_samples(lusc_bulk, r_hf)
lusc_bulk_low = filter_tx_samples(lusc_bulk, r_low)
lusc_bulk_mid = filter_tx_samples(lusc_bulk, r_mid)
lusc_bulk_high = filter_tx_samples(lusc_bulk, r_high)

## Make cache file listsings
wpure_lf_out = fs::path(outdir, 'wgbs_pure_lf.qs')
wpure_hf_out = fs::path(outdir, 'wgbs_pure_hf.qs')
wpure_low_out = fs::path(outdir, 'wgbs_pure_low.qs')
wpure_mid_out = fs::path(outdir, 'wgbs_pure_mid.qs')
wpure_high_out = fs::path(outdir, 'wgbs_pure_high.qs')
wbulk_lf_out = fs::path(outdir, 'wgbs_bulk_lf.qs')
wbulk_hf_out = fs::path(outdir, 'wgbs_bulk_hf.qs')
wbulk_low_out = fs::path(outdir, 'wgbs_bulk_low.qs')
wbulk_mid_out = fs::path(outdir, 'wgbs_bulk_mid.qs')
wbulk_high_out = fs::path(outdir, 'wgbs_bulk_high.qs')
luad_pure_lf_out = fs::path(outdir, 'luad_pure_lf.qs')
luad_pure_hf_out = fs::path(outdir, 'luad_pure_hf.qs')
luad_pure_low_out = fs::path(outdir, 'luad_pure_low.qs')
luad_pure_mid_out = fs::path(outdir, 'luad_pure_mid.qs')
luad_pure_high_out = fs::path(outdir, 'luad_pure_high.qs')
luad_bulk_lf_out = fs::path(outdir, 'luad_bulk_lf.qs')
luad_bulk_hf_out = fs::path(outdir, 'luad_bulk_hf.qs')
luad_bulk_low_out = fs::path(outdir, 'luad_bulk_low.qs')
luad_bulk_mid_out = fs::path(outdir, 'luad_bulk_mid.qs')
luad_bulk_high_out = fs::path(outdir, 'luad_bulk_high.qs')
lusc_pure_lf_out = fs::path(outdir, 'lusc_pure_lf.qs')
lusc_pure_hf_out = fs::path(outdir, 'lusc_pure_hf.qs')
lusc_pure_low_out = fs::path(outdir, 'lusc_pure_low.qs')
lusc_pure_mid_out = fs::path(outdir, 'lusc_pure_mid.qs')
lusc_pure_high_out = fs::path(outdir, 'lusc_pure_high.qs')
lusc_bulk_lf_out = fs::path(outdir, 'lusc_bulk_lf.qs')
lusc_bulk_hf_out = fs::path(outdir, 'lusc_bulk_hf.qs')
lusc_bulk_low_out = fs::path(outdir, 'lusc_bulk_low.qs')
lusc_bulk_mid_out = fs::path(outdir, 'lusc_bulk_mid.qs')
lusc_bulk_high_out = fs::path(outdir, 'lusc_bulk_high.qs')

## Run ORA for each file list, defining a different purity threshold
log_info("Running ORA for each purity threshold")
wgbs_lf_ora = pipe_ora_wgbs(wgbs_pure_lf, wgbs_bulk_lf, wpure_lf_out, wbulk_lf_out, fs::path(outdir, 'ORA_wgbs_lf.qs'))
wgbs_hf_ora = pipe_ora_wgbs(wgbs_pure_hf, wgbs_bulk_hf, wpure_hf_out, wbulk_hf_out, fs::path(outdir, 'ORA_wgbs_hf.qs'))
wgbs_low_ora = pipe_ora_wgbs(wgbs_pure_low, wgbs_bulk_low, wpure_low_out, wbulk_low_out, fs::path(outdir, 'ORA_wgbs_low.qs'))
wgbs_mid_ora = pipe_ora_wgbs(wgbs_pure_mid, wgbs_bulk_mid, wpure_mid_out, wbulk_mid_out, fs::path(outdir, 'ORA_wgbs_mid.qs'))
wgbs_high_ora = pipe_ora_wgbs(wgbs_pure_high, wgbs_bulk_high, wpure_high_out, wbulk_high_out, fs::path(outdir, 'ORA_wgbs_high.qs'))
luad_lf_ora = pipe_ora_rrbs(luad_pure_lf, luad_bulk_lf, luad_pure_lf_out, luad_bulk_lf_out, fs::path(outdir, 'ORA_luad_lf.qs'))
luad_hf_ora = pipe_ora_rrbs(luad_pure_hf, luad_bulk_hf, luad_pure_hf_out, luad_bulk_hf_out, fs::path(outdir, 'ORA_luad_hf.qs'))
luad_low_ora = pipe_ora_rrbs(luad_pure_low, luad_bulk_low, luad_pure_low_out, luad_bulk_low_out, fs::path(outdir, 'ORA_luad_low.qs'))
luad_mid_ora = pipe_ora_rrbs(luad_pure_mid, luad_bulk_mid, luad_pure_mid_out, luad_bulk_mid_out, fs::path(outdir, 'ORA_luad_mid.qs'))
luad_high_ora = pipe_ora_rrbs(luad_pure_high, luad_bulk_high, luad_pure_high_out, luad_bulk_high_out, fs::path(outdir, 'ORA_luad_high.qs'))
lusc_lf_ora = pipe_ora_rrbs(lusc_pure_lf, lusc_bulk_lf, lusc_pure_lf_out, lusc_bulk_lf_out, fs::path(outdir, 'ORA_lusc_lf.qs'))
lusc_hf_ora = pipe_ora_rrbs(lusc_pure_hf, lusc_bulk_hf, lusc_pure_hf_out, lusc_bulk_hf_out, fs::path(outdir, 'ORA_lusc_hf.qs'))
lusc_low_ora = pipe_ora_rrbs(lusc_pure_low, lusc_bulk_low, lusc_pure_low_out, lusc_bulk_low_out, fs::path(outdir, 'ORA_lusc_low.qs'))
lusc_mid_ora = pipe_ora_rrbs(lusc_pure_mid, lusc_bulk_mid, lusc_pure_mid_out, lusc_bulk_mid_out, fs::path(outdir, 'ORA_lusc_mid.qs'))
lusc_high_ora = pipe_ora_rrbs(lusc_pure_high, lusc_bulk_high, lusc_pure_high_out, lusc_bulk_high_out, fs::path(outdir, 'ORA_lusc_high.qs'))

# Generate complete ORA results table titled <ORA_RESF.tsv>
log_info("Generating complete ORA results table")
wgbs_lf = qs::qread(fs::path(outdir, 'ORA_wgbs_lf.qs'))
wgbs_hf = qs::qread(fs::path(outdir, 'ORA_wgbs_hf.qs'))
wgbs_low = qs::qread(fs::path(outdir, 'ORA_wgbs_low.qs'))
wgbs_mid = qs::qread(fs::path(outdir, 'ORA_wgbs_mid.qs'))
wgbs_high = qs::qread(fs::path(outdir, 'ORA_wgbs_high.qs'))
luad_lf = qs::qread(fs::path(outdir, 'ORA_luad_lf.qs'))
luad_hf = qs::qread(fs::path(outdir, 'ORA_luad_hf.qs'))
luad_low = qs::qread(fs::path(outdir, 'ORA_luad_low.qs'))
luad_mid = qs::qread(fs::path(outdir, 'ORA_luad_mid.qs'))
luad_high = qs::qread(fs::path(outdir, 'ORA_luad_high.qs'))
lusc_lf = qs::qread(fs::path(outdir, 'ORA_lusc_lf.qs'))
lusc_hf = qs::qread(fs::path(outdir, 'ORA_lusc_hf.qs'))
lusc_low = qs::qread(fs::path(outdir, 'ORA_lusc_low.qs'))
lusc_mid = qs::qread(fs::path(outdir, 'ORA_lusc_mid.qs'))
lusc_high = qs::qread(fs::path(outdir, 'ORA_lusc_high.qs'))
ora_res = data.table::rbindlist(
    list(
        wgbs_lf$ora %>% dplyr::mutate(source='mCRPC', purity='low'),
        wgbs_hf$ora %>% dplyr::mutate(source='mCRPC', purity='high'),
        wgbs_low$ora %>% dplyr::mutate(source='mCRPC', purity='low'),
        wgbs_mid$ora %>% dplyr::mutate(source='mCRPC', purity='mid'),
        wgbs_high$ora %>% dplyr::mutate(source='mCRPC', purity='high'),
        luad_lf$ora %>% dplyr::mutate(source='LUAD', purity='low'),
        luad_hf$ora %>% dplyr::mutate(source='LUAD', purity='high'),
        luad_low$ora %>% dplyr::mutate(source='LUAD', purity='low'),
        luad_mid$ora %>% dplyr::mutate(source='LUAD', purity='mid'),
        luad_high$ora %>% dplyr::mutate(source='LUAD', purity='high'),
        lusc_lf$ora %>% dplyr::mutate(source='LUSC', purity='low'),
        lusc_hf$ora %>% dplyr::mutate(source='LUSC', purity='high'),
        lusc_low$ora %>% dplyr::mutate(source='LUSC', purity='low'),
        lusc_mid$ora %>% dplyr::mutate(source='LUSC', purity='mid'),
        lusc_high$ora %>% dplyr::mutate(source='LUSC', purity='high')
    )
)
data.table::fwrite(ora_res, fs::path(outdir, 'ORA_RESF_PURITY.tsv'), sep='\t', quote=F)
