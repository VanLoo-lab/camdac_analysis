# Plot summary data at a single DMR
library(ggplot2)
library(data.table)
library(qs)
library(glue)
library(stringr)
library(ggcoverage)
library(GenomicRanges)
library(figpatch)
library(logger)
source("src/fora.R")
options(bitmapType = "cairo")

# Setup
ora_dir <- "data/pre/gsea_300"
outdir <- fs::dir_create("data/pre/gsea_300_top")
figdir = fs::dir_create(fs::path('figures/241006_fig_ora_top', "prom_figs"))
compfigdir=fs::dir_create(fs::path('figures/241006_fig_ora_top', "comp_figs"))

#  List which sample is in the `pure_lq` group for each dataset.
flist_lq <- list(
    "luad" = extract_sid_from_dmrlist(fs::path(ora_dir, "luad_pure_lq.qs")),
    "lusc" = extract_sid_from_dmrlist(fs::path(ora_dir, "lusc_pure_lq.qs")),
    "wgbs" = extract_sid_from_dmrlist(fs::path(ora_dir, "wgbs_pure_lq.qs"))
)

# Load methylation matrices
ro <- pipeline_extract_meth_matrix(flist_lq, outdir)
rmav <- get_group_averages(ro, outdir, cache='lq_group_averages.qs')

# Load recurrent DMRs
recc <- qs::qread(fs::path(outdir, "recc_lq.qs"))
#recc <- add_rmav_to_recc(recc, rmav)

# Load ORA results
ora_lq <- qs::qread(fs::path(outdir, "lq_top.qs"))

# Cache objects
ora_obs = list(recc=recc, rmav=rmav, ora_lq=ora_lq)
qs::qsave(ora_obs, fs::path(outdir, "ora_obs.qs"))

#  WGBS HYPERMETHYLATION ----
dset = 'wgbs'; dtype='hyper';
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal")
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# WGBS HYPOMETHYLATION ----
dset = 'wgbs'; dtype='hypo'; adtype=dtype
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal")
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::filter(dtype==adtype) %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# LUAD HYPERMETHYLATION ----
dset = 'luad'; dtype='hyper'; adtype=dtype
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal")
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::filter(dtype==adtype) %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# LUAD HYPOMETHYLATION ----
dset = 'luad'; dtype='hypo'; adtype=dtype
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal")
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::filter(dtype==adtype) %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# LUSC HYPERMETHYLATION ----
dset = 'lusc'; dtype='hyper'; adtype=dtype
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal")
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::filter(dtype==adtype) %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# LUSC HYPOMETHYLATION ----
dset = 'lusc'; dtype='hyper'; adtype=dtype
bulk=glue("{dset}_bulk"); pure=glue("{dset}_pure"); normal=glue("{dset}_normal") # 
sig_pway <- ora_lq$lq_pure[[dset]] %>% dplyr::filter(dtype==adtype) %>% dplyr::top_n(5, -pure) # Currently just top 5
gene_locs_all = recc[[dset]]$gene_data[[dtype]]$gene_loc # Gene locations for analysis
pipe_plot_spway_genes(sig_pway, recc, dset, dtype, rmav, figdir)

# Composites ----
# Use figpatch to combine the figures we want to show in the manuscript at specific genes
# Make one figure for each of the plot types (dot or line)

# Load data
fp_pcov = list(
    'luad'=c(
        'figures/241006_fig_ora_top/prom_figs/pcov_luad_hyper_KRAS.BREAST_UP.V1_UP_GFI1_6.png', # PMID 30606770
        'figures/241006_fig_ora_top/prom_figs/pcov_luad_hyper_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GREM1_4.png', # PMID 38306567
        'figures/241006_fig_ora_top/prom_figs/pcov_luad_hyper_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_LAMA2_10.png' # PMID 32792503
    ),
    'lusc'=c(
        'figures/241006_fig_ora_top/prom_figs/pcov_lusc_hyper_HALLMARK_PANCREAS_BETA_CELLS_PAX6_4.png', # PMID 26617874
        'figures/241006_fig_ora_top/prom_figs/pcov_lusc_hyper_JNK_DN.V1_UP_EFHD1_1.png',
        'figures/241006_fig_ora_top/prom_figs/pcov_lusc_hyper_JNK_DN.V1_UP_EFHD1_1.png'
    ),
    'mcrpc'=c(
        'figures/241006_fig_ora_top/prom_figs/pcov_wgbs_hyper_E2F3_UP.V1_UP_LEF1_1.png', ##
        'figures/241006_fig_ora_top/prom_figs/pcov_wgbs_hyper_CAMP_UP.V1_DN_DLG5_1.png',## Loss of Dlg5 expression promotes the migration and invasion of prstate cancer cells
        'figures/241006_fig_ora_top/prom_figs/pcov_wgbs_hyper_E2F3_UP.V1_UP_HEY2_5.png'
    )
)
fp_pdot = lapply(fp_pcov, function (x) str_replace(x, 'pcov_', 'pdot_'))

# Make figpatches for each dataset
make_fp <- function(x){
    patches = lapply(x, figpatch::fig)
    fpwrap = figpatch::fig_wrap(patches, ncol=1)
}
pcov_pat = lapply(fp_pcov, make_fp)
pcov_full = patchwork::wrap_plots(pcov_pat, nrow=1)
ggsave(fs::path(compfigdir, "pcov_full.png"), pcov_full, width=8, height=12)

pdot_pat = lapply(fp_pdot, make_fp)
pdot_full = patchwork::wrap_plots(pdot_pat, nrow=1)
ggsave(fs::path(compfigdir, "pdot_full.png"), pdot_full, width=8, height=12)
