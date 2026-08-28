# From ORA results, generate a data object containing:
# - [X] Top 10 pathways and genesets by significance in CAMDAC pure
# - [X] Top 10 genes by delta between CAMDAC pure and bulk
# - [X] IDs of samples in each of the ORA datasets (for lq only)
# - Genomic locations of each gene
# - Average DNAme matrix of bulk/pure/normal in a 1kb window around each gene

# This is so that we can highlight specific genes that would be discovered with CAMDAC alone

# Dependencies ----
library(fs)
library(dplyr)
library(stringr)
library(data.table)
library(logger)
library(ggplot2)
source('src/fora.R')
options(bitmapType = "cairo")

# Setup ----
ora_dir="data/pre/gsea_500"
outdir=fs::dir_create('data/pre/gsea_500/gsea_500_top')

## Top Pure and Top Pure Q-Delta ----
lq_ora = c(
    "luad"=fs::path(ora_dir, "ORA_luad_lq.qs"),
    "lusc"=fs::path(ora_dir, "ORA_lusc_lq.qs"),
    "wgbs"=fs::path(ora_dir, "ORA_wgbs_lq.qs")
)

## Adjust p-values if necessary
## update_pvals(lq_ora)

# Select most significant samples in all groups
lq_o10 = lapply(lq_ora, function(e){
    ora=parse_ora_f(e) %>%
        dplyr::group_by(dset, dtype, source) %>%
        dplyr::filter(padj <= 0.05)
})
lq_pure = lapply(lq_ora, function(e){
    # Select samples with shighest delta when N.S. in bulk (discovered by CAMDAC)
    ora = parse_ora_f(e) %>%
        dplyr::select(pathway, dset, dtype, padj, source, overlapGenes)
    odelta = ora %>%
        dplyr::select(pathway, dset, dtype, padj) %>% 
        tidyr::pivot_wider(names_from = dset, values_from = padj) %>% 
        dplyr::group_by(dtype) %>%
        dplyr::mutate(delta=abs(-log10(pure)--log10(bulk))) #%>% 
        #dplyr::filter(bulk > 0.05 & pure <= 0.05) 
    ora_merge = ora %>% dplyr::filter(dset=='pure') %>% dplyr::select(pathway, source, dtype, pure_edge=overlapGenes)
    left_join(odelta, ora_merge)
})
qs::qsave(list(lq_o10=lq_o10, lq_pure=lq_pure), fs::path(outdir, 'lq_top.qs'))

# Sample-level data ----
# List which sample is in the `pure_lq` group for each dataset.
flist_lq = list(
    "luad" = extract_sid_from_dmrlist(fs::path(ora_dir, 'luad_pure_lq.qs')),
    "lusc" = extract_sid_from_dmrlist(fs::path(ora_dir, 'lusc_pure_lq.qs')),
    "wgbs" = extract_sid_from_dmrlist(fs::path(ora_dir, 'wgbs_pure_lq.qs'))
)
# Summarise recurrent DMRs in each dataset
recc_lq <- lapply(
    c("luad", "lusc", "wgbs"),
    function(x){
        extract_recc_dmr(dset=x, ora_dir, outdir)
    }
)
names(recc_lq) <- c("luad", "lusc", "wgbs")
qs::qsave(recc_lq, fs::path(outdir, 'recc_lq.qs'))

