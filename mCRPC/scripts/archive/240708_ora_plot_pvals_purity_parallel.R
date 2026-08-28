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
figdir = fs::dir_create('figures/240628_ora_04')
suffix = c('lf', 'hf', 'lq', 'hq', 'low', 'mid', 'high')

for(si in suffix){
log_info("Generating ORA plots for {si}")
ora_files = fs::dir_ls(outdir, glob=paste0('*ORA_*', si, '.qs'))
ores = data.table::rbindlist(lapply(ora_files, parse_ora_f))

hmw = ores %>% 
    dplyr::mutate(lpadj2=-log10(padj2)) %>%
    dplyr::mutate(index=1:nrow(.)) %>%
    dplyr::select(pathway, dtype, hallmark, source, dset, lpadj2) %>%
    tidyr::pivot_wider(names_from=dset, values_from=lpadj2)
dout = fs::path(outdir, paste0('figd_ora_res_dp_', si, '.csv'))
data.table::fwrite(hmw, dout)

axscale = ifelse(hmw$dtype[[1]]=="hyper", 40, 3)
# Plot hypers
gphy = ggplot(hmw %>% dplyr::filter(dtype=='hyper'), aes(x=bulk, y=pure, fill=hallmark)) + 
    geom_vline(xintercept=-log10(0.05), color='red', alpha=0.8) +
    geom_hline(yintercept=-log10(0.05), color='red', alpha=0.8) +
    geom_point(shape=21, alpha=.6) + 
    facet_wrap(dtype~source, nrow=1) + geom_abline(intercept=0, slope=1, linetype='dashed') + theme_bw()

gpho = ggplot(hmw %>% dplyr::filter(dtype!='hyper'), aes(x=bulk, y=pure, fill=hallmark)) + 
    geom_vline(xintercept=-log10(0.05), color='red', alpha=0.8) +
    geom_hline(yintercept=-log10(0.05), color='red', alpha=0.8) +
    geom_point(shape=21, alpha=.6) + 
    facet_wrap(dtype~source, nrow=1) + geom_abline(intercept=0, slope=1, linetype='dashed') + theme_bw()

gp = patchwork::wrap_plots(gphy, gpho, ncol=1)
figout = fs::path(figdir, paste0('ora_res_dp_', si, '.png'))
ggsave(figout, gp, width=10, height=10, dpi=300)

}

