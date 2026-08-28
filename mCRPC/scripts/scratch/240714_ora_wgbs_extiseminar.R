# Aim: Plot p values using purity / parallel
library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(fs)
library(qs)
library(logger)
library(patchwork)
options(bitmapType='cairo')

# Functions ----

parse_ora_f <- function(f){
    e = fs::path_file(f) %>% stringr::str_remove('.qs')
    e = stringr::str_split(e, '_') %>% unlist
    dt = qs::qread(f)$ora
    dt$source=e[[2]]
    dt$purityf=e[[3]]
    return(dt)
}

# Main ----
outdir = fs::dir_create("data/pre/gsea_04")
figdir = fs::dir_create('figures/240714_ora_exitsem')
suffix = c('lf', 'hf', 'lq', 'hq', 'low', 'mid', 'high')
suffix = c('lq', 'hq')

for(si in suffix){
log_info("Generating ORA plots for {si}")
ora_files = fs::dir_ls(outdir, glob=paste0('*ORA_*', si, '.qs'))
ores = data.table::rbindlist(lapply(ora_files, parse_ora_f))
ores = dplyr::filter(ores, source=="wgbs")

hmw = ores %>% 
    dplyr::mutate(lpadj2=-log10(padj2)) %>%
    dplyr::mutate(index=1:nrow(.)) %>%
    dplyr::select(pathway, dtype, hallmark, source, dset, lpadj2) %>%
    tidyr::pivot_wider(names_from=dset, values_from=lpadj2) %>%
    dplyr::mutate(dtype= ifelse(dtype=="hyper", "Hypermethylation", "Hypomethylation"))

# Plot hypers
gphy = ggplot(hmw , aes(x=bulk, y=pure, fill=hallmark)) + 
    geom_point(shape=21, alpha=.5) + 
    facet_wrap(~dtype, ncol=1, scales="free") + geom_abline(intercept=0, slope=1, linetype='dashed') + theme_bw() +
    theme(strip.background = element_blank()) +
    labs(x="Bulk tumour\nORA q-value (-log10)", y="Pure tumour\nover-representation q-value (-log10)", fill="MSigDb Collection") +
    scale_fill_manual(labels=c("hall"="Hallmark Pathways", "o3ca"="Curated Cancer\nCell Atlas", "onc"="Oncogenic\nSignature Genes"), values=c("black", "orange", "purple"))

figout = fs::path(figdir, paste0('ora_res_dp_', si, '.png'))
ggsave(figout, gphy, width=5, height=8, dpi=300)
figout = fs::path(figdir, paste0('ora_res_dp_', si, '.pdf'))
ggsave(figout, gphy, width=5, height=8, dpi=300)

}

