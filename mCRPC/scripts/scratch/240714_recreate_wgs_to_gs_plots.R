# Recreate agg_pipe plots
library(qs)
library(ggplot2)
options(bitmapType = 'cairo')


plot_intersects <- function(e){
    e = dplyr::arrange(e, desc(perc_overlap))
    e$samples = factor(seq(nrow(e)))
    x= tidyr::pivot_longer(e, c("perc_overlap", "perc_overlap_1"), names_to="overlap_type", values_to="perc_overlap")
    x$overlap_type = ifelse(x$overlap_type == "perc_overlap", "exact", "within_1")
    x = dplyr::filter(x, overlap_type=="exact")
    xax=paste0(x$name_a[1], "", x$method_a[1])
    yax=paste0(x$name_b[1], "", x$method_b[1])
    ggplot(x, aes(x=samples, y=perc_overlap, fill=overlap_type)) +
        geom_point(shape=23, color='black') +
        scale_fill_manual(values=c("darkgreen", "orange")) +
        labs(x="Sample", y="Exact match (% genome)", title="WGBS vs WGS Copy number overlap") +
        theme_classic() +
        scale_y_continuous(limits=c(0,100)) +
        theme(axis.text.x=element_blank()) +
        # place legend within plot at bottom right
        theme(legend.position="none")

}

plot_ploidies <- function(x){
    xax=paste0(x$name_a[1], "", x$method_a[1])
    yax=paste0(x$name_b[1], "", x$method_b[1])
    x$ploidy_a = ifelse(x$ploidy_a > 6, 6, x$ploidy_a)
    x$ploidy_b = ifelse(x$ploidy_b > 6, 6, x$ploidy_b)
    # Set last label on X and Y axis to '6+'
    xax_labels = seq(0,6,1)
    xax_labels[length(xax_labels)] = "6+"
    yax_labels = seq(0,6,1)
    yax_labels[length(yax_labels)] = "6+"
    ggplot(x, aes(x=ploidy_a, y=ploidy_b)) +
    geom_abline(intercept=0, slope=1, color='grey', linetype='dashed') +
    geom_abline(intercept=2, slope=1, color='grey', linetype='dashed') +
    geom_abline(intercept=-2, slope=1, color='grey', linetype='dashed') +
    geom_point(shape=23, color='black', fill='skyblue') +
    scale_x_continuous(limits=c(0,6), breaks=seq(0,6,1), labels=xax_labels) +
    scale_y_continuous(limits=c(0,6), breaks=seq(0,6,1), labels=yax_labels) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x=xax, y=yax, title="Ploidy") +
    theme_classic()
}

plot_purities <- function(x){
    xax=paste0(x$name_a[1], "", x$method_a[1])
    yax=paste0(x$name_b[1], "", x$method_b[1])
    ggplot(x, aes(x=purity_a, y=purity_b)) + 
    geom_abline(intercept=0, slope=1, color='grey', linetype='dashed') +
    geom_point(shape=23, color='black', fill='purple') +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x=xax, y=(yax), title="Purity") +
    theme_classic()
}


ppplot <- function(data, outfile){
    p_purity = plot_purities(data)
    p_ploidy = plot_ploidies(data)
    p_intersect = plot_intersects(data)
    p_wrap = patchwork::wrap_plots(p_purity, p_ploidy, p_intersect, nrow=1)
    ggsave(outfile, p_wrap, width=12, height=4)
}

resdir = fs::dir_create('figures/240714_wgs_to_gs')
to = qs::qread("data/archive/2023_data/230916_wgbs_ascat/intersects/bsto_wgs_bb.qs")$res
to$method_a=""; to$method_b=""; to$name_a="WGS Tumour-Normal"; to$name_b="WGBS Tumour Only"
ppplot(to, fs::path(resdir, "bsto_wgs_bb.png"))
ppplot(to, fs::path(resdir, "bsto_wgs_bb.pdf"))
gs = qs::qread("data/archive/2023_data/230916_wgbs_ascat/intersects/bsgs_wgs_bbn.qs")$res
gs$method_a=""; gs$method_b=""; gs$name_a="WGS Tumour-Normal"; gs$name_b="WGBS Tumour Only with known\nGermline Het. SNP position"
ppplot(gs, fs::path(resdir, "bsgs_wgs_bbn.png"))
ppplot(gs, fs::path(resdir, "bsgs_wgs_bbn.pdf"))
