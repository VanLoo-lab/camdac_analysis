# Plot over-enrichment analysis results
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
options(bitmapType='cairo')

clean_pathways <- function(x) stringr::str_remove(x, '^GAVISH_3CA_MALIGNANT_METAPROGRAM_') %>% gsub('METAPROGRAM_', '', .) %>% gsub('HALLMARK_', '', .)

outdir = fs::dir_create('figures/240416_ora')
ora = 'data/pre/gsea/ora_res.tsv'
gora = fread(ora) %>% dplyr::mutate(perc_ol=overlap*100/size, dset=str_to_title(dset), dtype = str_to_title(dtype), pathway=clean_pathways(pathway))

# Convert padj2 field to stars for signifiance
gora$padj_stars = cut(gora$padj2, breaks=c(0, 0.001, 0.01, 0.05, 1), labels=c('***', '**', '*', ''))

hora = gora %>% dplyr::filter(hallmark=='hall')
aora = gora %>% dplyr::filter(hallmark=='o3ca')
cora = gora %>% dplyr::filter(hallmark=='onc')
cora_nb = cora %>% dplyr::filter(hallmark=='onc' & pathway %in% (dplyr::filter(cora, dset=='Pure') %>% pull(pathway)))

# Capitalise the dset field in gora

plot_ora_bars <- function(x){
ggplot(x, aes(x=perc_ol, y=pathway, fill=dset)) +
    geom_col(position = position_dodge()) + 
    facet_grid(dtype ~ source, scales='free_y', space='free') +
    scale_x_continuous(limits=c(0,100)) +
    # Add geom_text labl with pvalue starts above the bar
    geom_text(aes(label=padj_stars, group=dset), position=position_dodge(width=0.75), size=3, hjust=-0.4) +
    theme_bw() +
    labs(x="Genes (%)", y="", fill="Source") +
    scale_fill_manual(values=c('skyblue', 'purple')) +
    theme(legend.position = 'bottom', strip.background = element_blank())
}

# Generate plots
p_hora = plot_ora_bars(hora) + labs(title="Hallmark (MSigDb:H)")
p_aora = plot_ora_bars(aora) + labs(title="GAVISH Malignant Metaprogram (MSigDb:C4:3CA)")
p_cora = plot_ora_bars(cora) + labs(title="Oncogenic (MSigDb:C6)")
p_cora_nb = plot_ora_bars(cora_nb) + labs(title="Oncogenic (MSigDb:C6) Pure-Filtered")

ggsave(fs::path(outdir, 'hallmark_ora.png'), plot_ora_bars(hora), width=8, height=12)
ggsave(fs::path(outdir, 'aora_ora.png'), plot_ora_bars(aora), width=8, height=12)
ggsave(fs::path(outdir, 'cora_ora.png'), plot_ora_bars(cora), width=8, height=12) # Many LUSC signals resulting from bulk only
ggsave(fs::path(outdir, 'cora_nb_ora.png'), plot_ora_bars(cora_nb), width=8, height=12)
