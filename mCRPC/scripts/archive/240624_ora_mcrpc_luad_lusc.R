
# Aim: Over-enrichment analysis on bulk and pure promoter DMRs for mCRPC, LUAD, and LUSC samples
# ORA analysis tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html

# Setup -----
source('src/fora.R')
library_ora()
library(logger)
library(ggplot2)

# Files and directories
outdir = fs::dir_create("data/pre/gsea")
wgbs_root = "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/dmrs/"
rrbs_root = "/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMRs/"

# Run WGBS DMR analysis
wpure_out = fs::path(outdir, 'wgbs_pure_DMR_list.qs')
wbulk_out = fs::path(outdir, 'wgbs_bulk_DMR_list.qs')
wgbs_ora_out = fs::path(outdir, 'wgbs_prom_DMR_ora.qs')
wgbs_pure = fs::dir_ls(wgbs_root, recurse=T, glob="*pure_dmr.fst")
wgbs_bulk= fs::dir_ls(wgbs_root, recurse=T, glob="*bulk_dmr.fst")
wgbs_ora = pipe_ora_wgbs(wgbs_pure, wgbs_bulk, wpure_out, wbulk_out, wgbs_ora_out)

# Run RRBS DMR analysis for LUAD
luad_pure_out = fs::path(outdir, 'rrbs_luad_pure_DMR_list.qs')
luad_bulk_out = fs::path(outdir, 'rrbs_luad_bulk_DMR_list.qs')
luad_ora_out = fs::path(outdir, 'rrbs_luad_prom_DMR_ora.qs')
luad_pure = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUAD")
luad_bulk = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUAD")
luad_ora = pipe_ora_rrbs(luad_pure, luad_bulk, luad_pure_out, luad_bulk_out, luad_ora_out)

# Run RRBS DMR analysis for LUSC
lusc_pure_out = fs::path(outdir, 'rrbs_lusc_pure_DMR_list.qs')
lusc_bulk_out = fs::path(outdir, 'rrbs_lusc_bulk_DMR_list.qs')
lusc_ora_out = fs::path(outdir, 'rrbs_lusc_prom_DMR_ora.qs')
lusc_pure = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*purified_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUSC")
lusc_bulk = fs::dir_ls(rrbs_root, recurse=T, glob="*re_segmented*bulk_tumour_DMR_calls.RData") %>% keep_rrbs_subtype(subt="LUSC")
lusc_ora = pipe_ora_rrbs(lusc_pure, lusc_bulk, lusc_pure_out, lusc_bulk_out, lusc_ora_out)

# Generate dataset for plotting, should contain pathway results for all three analyses
wgbs_ora = qs::qread(wgbs_ora_out)
luad_ora = qs::qread(luad_ora_out)
lusc_ora = qs::qread(lusc_ora_out)
ora_res = data.table::rbindlist(
    list(
        wgbs_ora$ora %>% dplyr::mutate(source='mCRPC'),
        luad_ora$ora %>% dplyr::mutate(source='LUAD'),
        lusc_ora$ora %>% dplyr::mutate(source='LUSC')
    )
)
data.table::fwrite(ora_res, fs::path(outdir, 'ora_res.tsv'), sep='\t', quote=F)
