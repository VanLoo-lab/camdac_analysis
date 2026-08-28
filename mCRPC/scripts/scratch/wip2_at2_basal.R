library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(patchwork)
figdir = fs::dir_create('figures/240703_at2_basal')
options(bitmapType='cairo')

# Load dataset and perform cleanup operations
res = qs::qread('data/pre/at2_basal/at2_basal_regions_v2.qs')
res = res %>% dplyr::mutate(gene=ifelse(rid==12, 'TP63_RNA', gene)) 
res = res %>% dplyr::filter(dset=='genes')
#res = res %>% dplyr::filter(gene %in% c("SFTPA1", "SFTPB", "SFTPC", "KRT5", "TP63", "GJB5"))

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')

cpure = ggplot(res, aes(x=gene, y=mt_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="CAMDAC gene body methylation", x="") + 
    scale_y_continuous(limits=c(0,1)) 

cbulk = ggplot(res, aes(x=gene, y=mb_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="Bulk gene body methylation", x="") +
    scale_y_continuous(limits=c(0,1)) 

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')
ggsave('figures/240703_at2_basal/fig_at2_basal_genes.png', pfig, width=12, height=6)

# Load dataset and perform cleanup operations
res = qs::qread('data/pre/at2_basal/at2_basal_regions_v2.qs')
res = res %>% dplyr::mutate(gene=ifelse(rid==12, 'TP63_RNA', gene)) 
res = res %>% dplyr::filter(dset=='promoters')
#res = res %>% dplyr::filter(gene %in% c("SFTPA1", "SFTPB", "SFTPC", "KRT5", "TP63", "GJB5"))

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')

cpure = ggplot(res, aes(x=gene, y=mt_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="CAMDAC gene body methylation", x="") + 
    scale_y_continuous(limits=c(0,1)) 

cbulk = ggplot(res, aes(x=gene, y=mb_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="Bulk gene body methylation", x="") +
    scale_y_continuous(limits=c(0,1)) 

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')
ggsave('figures/240703_at2_basal/fig_at2_basal_promoters.png', pfig, width=8, height=6)

# Load dataset and perform cleanup operations
res = qs::qread('data/pre/at2_basal/at2_basal_regions_v2.qs')
res = res %>% dplyr::mutate(gene=ifelse(rid==12, 'TP63_RNA', gene)) 
res = res %>% dplyr::filter(dset=='enhancers')
#res = res %>% dplyr::filter(gene %in% c("SFTPA1", "SFTPB", "SFTPC", "KRT5", "TP63", "GJB5"))

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')

cpure = ggplot(res, aes(x=gene, y=mt_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="CAMDAC gene body methylation", x="") + 
    scale_y_continuous(limits=c(0,1)) 

cbulk = ggplot(res, aes(x=gene, y=mb_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="Bulk gene body methylation", x="") +
    scale_y_continuous(limits=c(0,1)) 

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')
ggsave('figures/240703_at2_basal/fig_at2_basal_enhancers.png', pfig, width=8, height=6)

