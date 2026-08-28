# Export ORA results for supplementary table 4
library(dplyr)
library(data.table)
source("src/fora.R")

# Ora
ora_dir <- "data/pre/gsea_500"

# Setup function to extract
ora_qread <- function(x, label, dq){
    ora = qs::qread(x)$ora
    ora = ora %>% dplyr::filter(hallmark=='onc') %>% 
        dplyr::mutate(dsource=label, dq=dq)
    ora %>% dplyr::select(
        pathway, padj, overlap, size, dtype, dset, hallmark, dsource, dq
    )
}

# Load file list
lq <- list(
    "wgbs" = ora_qread(fs::path(ora_dir, "ORA_wgbs_lq.qs"), "mcrpc", 'lq'),
    "lusc" = ora_qread(fs::path(ora_dir, "ORA_lusc_lq.qs"), "lusc", 'lq'),
    "luad" = ora_qread(fs::path(ora_dir, "ORA_luad_lq.qs"), "luad", 'lq')
) %>% data.table::rbindlist()
hq <- list(
    "wgbs" = ora_qread(fs::path(ora_dir, "ORA_wgbs_hq.qs"), "mcrpc", 'hq'),
    "lusc" = ora_qread(fs::path(ora_dir, "ORA_lusc_hq.qs"), "lusc", 'hq'),
    "luad" = ora_qread(fs::path(ora_dir, "ORA_luad_hq.qs"), "luad", 'hq')
) %>% data.table::rbindlist()

# Combine
xport = rbind(lq, hq)

xbp = xport %>% dplyr::select(-overlap, -size) %>%
    tidyr::pivot_wider(names_from=c(dset, dsource), values_from=padj) %>%
    dplyr::arrange(dtype, dq)

# Only show genes where one dataset is significant
xany_signif = xbp %>% dplyr::filter(
    bulk_mcrpc <= 0.05 | pure_mcrpc <= 0.05 | bulk_lusc <= 0.05 | pure_lusc <= 0.05 | bulk_luad <= 0.05 | pure_luad <= 0.05
) 
data.table::fwrite(xany_signif, 'data/out/st4_ora_export.tsv', sep='\t', quote=F, row.names=F)
# Report hypo: xport %>% dplyr::filter(dtype=='hypo' & padj<=0.05)

# What percentage of significant pathways are exclusive or shared?
xf = xany_signif %>% dplyr::filter(dtype=='hyper' & dq=='lq')
xc = xf %>% transmute(
    isbulk_mcrpc = bulk_mcrpc <= 0.05 & pure_mcrpc > 0.05,
    ispure_mcrpc = bulk_mcrpc > 0.05 & pure_mcrpc <= 0.05,
    isbulk_lusc = bulk_lusc <= 0.05 & pure_lusc > 0.05,
    ispure_lusc = bulk_lusc > 0.05 & pure_lusc <= 0.05,
    isbulk_luad = bulk_luad <= 0.05 & pure_luad > 0.05,
    ispure_luad = bulk_luad > 0.05 & pure_luad <= 0.05
)
xc_atleast_one = dplyr::mutate(xc,
    pure_hit = dplyr::case_when(
        ispure_mcrpc & !isbulk_mcrpc ~ T,
        ispure_lusc & !isbulk_lusc ~ T,
        ispure_luad & !isbulk_luad ~ T,
        TRUE ~ F
    ),
    bulk_hit = dplyr::case_when(
        isbulk_mcrpc & !ispure_mcrpc ~ T,
        isbulk_lusc & !ispure_lusc ~ T,
        isbulk_luad & !ispure_luad ~ T,
        TRUE ~ F
    ))
xc_atleast_one %>% colSums() / nrow(xc) * 100
# isbulk_mcrpc ispure_mcrpc  isbulk_lusc  ispure_lusc  isbulk_luad  ispure_luad 
#     2.500000    43.333333     3.333333    16.666667     1.666667    25.833333 
#     pure_hit     bulk_hit 
#    68.333333     6.666667 

# Next, give most striking examples
xfq = dplyr::filter(xf, dtype=='hyper' & dq=='lq')
dplyr::arrange(xf %>% filter(pure_mcrpc<=0.05), pure_mcrpc-bulk_mcrpc) %>% head
dplyr::arrange(xf  %>% filter(pure_lusc<=0.05), pure_lusc-bulk_lusc) %>% head
dplyr::arrange(xf  %>% filter(pure_luad<=0.05), pure_luad-bulk_luad) %>% head %>% data.table
