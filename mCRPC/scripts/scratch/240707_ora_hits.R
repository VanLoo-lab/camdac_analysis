# Aim: Are any of the ORA results biologically meaningful?
library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(fs)
library(qs)
library(logger)
library(gt)
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

ora_files = fs::dir_ls(outdir, glob=paste0('*ORA_*', 'lq', '.qs'))
ores = data.table::rbindlist(lapply(ora_files, parse_ora_f)) %>% 
    dplyr::mutate(lpadj2=-log10(padj2)) 
vores = data.table::rbindlist(lapply(ora_files, parse_ora_f)) %>% 
    dplyr::mutate(lpadj2=-log10(padj2)) %>%
    dplyr::mutate(index=1:nrow(.)) %>%
    dplyr::select(pathway, dtype, hallmark, source, dset, lpadj2) %>%
    tidyr::pivot_wider(names_from=dset, values_from=lpadj2)

vla = vores %>% dplyr::filter(source=='luad') %>% dplyr::arrange(-pure)
vlu = vores %>% dplyr::filter(source=='lusc') %>% dplyr::arrange(-pure)
vwg = vores %>% dplyr::filter(source=='wgbs') %>% dplyr::arrange(-pure)

# Plot table of top 10 per gene list for lq only
dplyr::filter(vwg, dtype=="hyper") %>% dplyr::arrange(-pure) %>% dplyr::filter(hallmark=="o3ca") %>% head(10)
dplyr::filter(vwg, dtype=="hyper") %>% dplyr::arrange(-pure) %>%dplyr::filter(hallmark=="o3ca" & bulk <= -log10(0.01)) %>% head(10)
dplyr::filter(vwg, dtype=="hypo") %>% dplyr::arrange(-pure) %>% dplyr::filter(hallmark=="o3ca") %>% head(10)
dplyr::filter(vwg, dtype=="hypo") %>% dplyr::arrange(-pure) %>%dplyr::filter(hallmark=="o3ca" & bulk <= -log10(0.01)) %>% head(10)

# Investigations
aj = dplyr::filter(ores, pathway=="HALLMARK_APICAL_JUNCTION" & source=='luad')
intersect(aj$overlapGenes[[1]], aj$overlapGenes[[2]])

head(aj)

