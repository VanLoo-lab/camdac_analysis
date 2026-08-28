library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsignif)
library(patchwork)
figdir = fs::dir_create('figures/240417_at2_basal')
options(bitmapType='cairo')

# Load dataset and perform cleanup operations
res = qs::qread('data/pre/at2_basal/at2_basal_regions.qs')
res = res %>% dplyr::mutate(gene=ifelse(rid==12, 'TP63_RNA', gene)) 
res = res %>% dplyr::filter(dset=='genes')
res = res %>% dplyr::filter(gene %in% c("SFTPA1", "SFTPB", "SFTPC", "KRT5", "TP63", "GJB5"))

# Calculate significance scores and padj for each gene between luad and lusc for bulk and pure
bsignif = res %>% group_by(gene,cell) %>% dplyr::summarize(
    pval=stats::wilcox.test(mb_av ~ subtype, alternative='two.sided')$p.value
) %>% dplyr::ungroup() %>% dplyr::mutate(dataset='bulk')
psignif = res %>% group_by(gene,cell) %>% dplyr::summarize(
    pval=stats::wilcox.test(mt_av ~ subtype, alternative='two.sided')$p.value
) %>% dplyr::ungroup() %>% dplyr::mutate(dataset='pure')
signif = rbind(bsignif, psignif) %>% dplyr::mutate(padj=p.adjust(pval, method='BH')) %>%
    dplyr::mutate(padj_stars = dplyr::case_when(
        padj < 0.001 ~ '***',
        padj < 0.01 ~ '**',
        padj < 0.05 ~ '*',
        TRUE ~ ''
    ))

cpure = ggplot(res, aes(x=gene, y=mt_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="CAMDAC gene body methylation", x="") + 
    scale_y_continuous(limits=c(0,1)) +
    # Add signif
    geom_signif(
        data = signif %>% dplyr::filter(dataset=="pure"),
        mapping=aes(xmin=gene, xmax=gene, annotations=sprintf("%.2e",padj), y_position=0.95),
        textsize=3, vjust=-.2, manual=TRUE, inherit.aes=FALSE, tip_length=c(0,0)
    )

cbulk = ggplot(res, aes(x=gene, y=mb_av, color=subtype)) + geom_boxplot(outlier.shape=NA) + 
    geom_jitter(shape=21, alpha=0.5, position=position_jitterdodge()) +
    facet_grid(~cell, scales='free_x', space='free') +
    theme_bw() + theme(strip.background=element_rect(fill='white')) +
    scale_color_manual(values=c("purple", "orange")) +
    labs(color='NSCLC', y="Bulk gene body methylation", x="") +
    scale_y_continuous(limits=c(0,1)) +
    # Add signif
    geom_signif(
        data = signif %>% dplyr::filter(dataset=="bulk"),
        mapping=aes(xmin=gene, xmax=gene, annotations=sprintf("%.2e",padj), y_position=0.95),
        textsize=3, vjust=-.2, manual=TRUE, inherit.aes=FALSE, tip_length=c(0,0)) 

pfig = wrap_plots(list(cbulk, cpure), ncol=1) + plot_layout(guides='collect')
ggsave('figures/240417_at2_basal/fig_at2_basal.png', pfig, width=8, height=6)
pfig
