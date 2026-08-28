library(data.table)
library(dplyr)
library(ggplot2)
library(logger)
library(patchwork)
library(ComplexHeatmap)
library(GenomicRanges)
library(scales)#
library(tidyr)
#library(gggenes)
options(bitmapType='cairo')

count_hmrs <- function(dfile){
    data = qs::qread(dfile)
    ccount = data$rHMRc[, .(ncg_meth_inhmrs=sum(ncg_meth, na.rm=T)), by='variable']
    names(ccount) = c("sample_id", "ncg_meth_inhmrs")
    return(ccount)
}

load_rhmrs <- function(dfile){

    # Data for 1 HMR cohort (bulk or pure)
    log_info("Reading HMR data")
    data = qs::qread(dfile)

    log_info("Pivoting HMR Data")
    # Pivot the data using dcast/melt
    rHMRwide <- data.table::dcast(
        data$rHMRc,
        rHMRix ~ variable,
        value.var = "mean_meth"
    )
    # Add nsamples with methylation data as anno to rHMRv
    rhmrcounts = data.table(
        rHMRix = rHMRwide[[1]],
        nsamp = rowSums(!is.na(rHMRwide[, -1]))
    )
    data$rHMRv = merge(data$rHMRv, rhmrcounts, by='rHMRix', all.x=T)
    data$rHMRwide = rHMRwide

    return(data)
}

get_top_rhmrs <- function(x, N){
    # Get rHMRs with data in atleast 25% of all samples
    nsamps = ncol(x$rHMRwide[, -1])
    nsamp_filt = round(nsamps * 0.25, 0)

    # Select rHMRs with coverage in the minimum number of samples
    rvd = x$rHMRv[nsamp >= nsamp_filt,]
    rvd = rvd[order(-rvd$rvar),]

    # Select top N rHMRs from the rvar rank
    topix = rvd[1:N]$rHMRix
    rtops = x$rHMRwide[rHMRix %in% topix,]
    rtops = rtops[,-1]

    # Create an output object
    obj = list(
        rtops_ix = topix,
        rtops = rtops,
        rvd = rvd
    )
    return(obj)
}

build_cg_map <- function(dfile, meta, cgi){

    # Load data
    x = qs::qread(dfile)

    # Count CpGs methylated in rHMRs
    ccount = count_hmrs(dfile)

    # Count CpGs by Cgi Annoataion
    hmr_cgi = count_hmr_cgi(x, cgi)

    # Next, generate the HMR boxplots
    bdata = merge(x$hmatc, meta, by='sample_id', all.x=T)
    bdata = merge(bdata, x$hmr_ps, by='sample_id', all.x=T)
    bdata = merge(bdata, ccount, by='sample_id', all.x=T)
   
    return(list(bdata=bdata, hmr_cgi=hmr_cgi))
}

count_hmr_cgi <- function(x, cgi){
    # Count CpGs in CGI regions
    cggr = GRanges(
        seqnames=cgi$chrom,
        ranges=IRanges(start=cgi$start, end=cgi$end),
        cgi_annot=cgi$cgi_annot
    )
    seqlevelsStyle(cggr) = 'UCSC'
    cgol = findOverlaps(x$rHMRs, cggr)
    xrd = data.table(
        rHMRix = x$rHMRs[queryHits(cgol)]$rHMRix,
        cgi_annot = cggr[subjectHits(cgol)]$cgi_annot
    )
    xrc = xrd[, .(p_osea=sum(cgi_annot=='opensea'), p_shelf=sum(cgi_annot=='shelf'), p_shore=sum(cgi_annot=='shore'), p_island=sum(cgi_annot=='islands'), ncg=.N), by='rHMRix']
    xrc$frac_osea = xrc$p_osea / xrc$ncg
    return(xrc)

}

count_bp_cgi <- function(x, cgi){
    hmr_mat = x$hmr_mat
    # Convert x to an object for the count
    x = list(rHMRs=x$gpos)
    x$rHMRs$rHMRix = seq(length(x$rHMRs))

    # Count CpGs in CGI regions
    ccount = count_hmr_cgi(x, cgi)
    ccount$gpos_ix = ccount$rHMRix
    ccount$rHMRix=NULL
    ccount$is_osea = ccount$frac_osea > 0.8

    # For each sample, count the number of HMRs
    hmr_mat$gpos_ix = seq(nrow(hmr_mat))
    hmr_mat = merge(hmr_mat, ccount[, .(is_osea, gpos_ix)], by='gpos_ix', all.x=T)
    hmr_mat[is.na(is_osea), is_osea:=T]
    n_osea = hmr_mat[is_osea==T,] %>% dplyr::select(-gpos_ix, -is_osea) %>% colSums
    nsamp = hmr_mat %>% dplyr::select(-gpos_ix, -is_osea) %>% colSums

    res = data.table(
        samp = names(nsamp),
        n_osea_100bp = n_osea,
        n_samp_100bp = nsamp
    )
    res$prop_osea = res$n_osea_100bp / res$n_samp_100bp

    return(res)

}

find_rhmr_sets <- function(brhmr, prhmr, brhtop, prhtop){
    # Generate a vector for all rHMRs in bulk and pure. Determine if they are in the other set by overlap and annotation. Annotate if they are amongst the most variable rHMrs (brhtop)

    bulkhrs = brhmr$rHMRs
    purehrs = prhmr$rHMRs

    # Annotate rHMRs by overlaps
    ols = findOverlaps(bulkhrs, purehrs)
    bol = mcols(bulkhrs) %>% as.data.frame() %>% data.table() %>% dplyr::mutate(dataset='bulk')
    suppressWarnings(bol[queryHits(ols), in_other := TRUE])
    pol = mcols(purehrs) %>% as.data.frame() %>% data.table() %>% dplyr::mutate(dataset='pure')
    suppressWarnings(pol[subjectHits(ols), in_other := TRUE])
    
    # Add annotation for variable rHMRs
    bol = dplyr::mutate(bol, in_top = rHMRix %in% brhtop$rtops_ix, in_other_top = rHMRix %in% prhtop$rtops_ix)
    pol = dplyr::mutate(pol, in_top = rHMRix %in% prhtop$rtops_ix, in_other_top = rHMRix %in% brhtop$rtops_ix)

    # Combine datasets
    rhmrd = rbind(bol, pol)
    rhmrd[is.na(in_other), in_other := FALSE]
    return(rhmrd)
}

plot_rhmr_heatmap <- function(tops_obj, hmr_sets, meta, hmr_cgi, figout){
    rtops = tops_obj$rtops
    qcpass = dplyr::filter(meta, tPASS==TRUE) %>% dplyr::pull(sample_id)
    rtops = data.table(rtops)[, ..qcpass]
    rtops_mat = as.matrix(rtops)
    rtops_mat[is.na(rtops_mat)] = 0.5

    # Setup annotations based on the metadata table
    rtmeta = meta[match(names(rtops), meta$sample_id),]
    rtmeta_s = dplyr::select(rtmeta, -sample_id, -tPASS, -cimp)

    # Setup row splits for cgi
    rcgi = data.table(rHMRix = tops_obj$rtops_ix)
    rcgi = merge(rcgi, hmr_cgi, by='rHMRix', all.x=T)
    osea = dplyr::case_when(
        is.na(rcgi$frac_osea) ~ "unknown",
        rcgi$frac_osea > 0.25 ~ 'opensea',
        rcgi$p_island > rcgi$p_shore + rcgi$p_shelf ~ 'island',
        rcgi$p_island > 0 ~ 'island',
        TRUE ~ 'shore/shelf'
    )

    # Setup row annotations for HMR sets
    rsets = data.table(rHMRix = tops_obj$rtops_ix) 
    rsets = merge(rsets, hmr_sets, by='rHMRix', all.x=T)
    rsets = dplyr::mutate(rsets, setclass=dplyr::case_when(
        in_top==TRUE & in_other_top==TRUE ~ "rHMR (both)",
        in_top==TRUE & in_other_top==FALSE ~ "rHMR (bulk)",
        in_top==FALSE & in_other_top==TRUE ~ "rHMR (pure)",
        TRUE ~ "rHMR (unknown)"
    ) )
    rsets = rsets$setclass


    # Create complex heatmap from rtops using data$hmrpat as annotation
    log_info("Building heatmap")
    hm = Heatmap(
        rtops_mat,
        name = "Methylation",
        column_title = "Samples",
        row_title = "rHMRs",
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(
            df = rtmeta_s,
            annotation_legend_param = list(
                purity=list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 13)),
                ploidy=list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 13)),
                stype=list(labels_gp =  gpar(fontsize = 13), title_gp = gpar(fontsize = 13)),
                tissue=list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 13)),
                nrpcc=list(labels_gp =  gpar(fontsize = 13), title_gp = gpar(fontsize = 13))
            )
        ),
        right_annotation = rowAnnotation(
            class=rsets,
            cgi=osea,
            annotation_legend_param = list(
                class=list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize=13)),
                cgi=  list(labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize=13))
            )
        ),
        row_split = osea,
        heatmap_legend_param = list(
            title = "M", at = c(0, 0.5, 1),
            labels_gp = gpar(fontsize = 13),
            title_gp = gpar(fontsize = 13)
            )
    )

    # Generate heatmap
    png(figout, width=12, height=12, units='in', res=300)
    draw(hm)
    dev.off()

    # Return figure and obs
    obj = list(
        hm = hm,
        rtmeta_s = rtmeta_s,
        rtops_mat = rtops_mat
    )
    return(hm)
}

plot_hcounts <- function(brcg, prcg) {
    hcounts = rbind(
        dplyr::mutate(brcg$bdata, dataset='bulk'),
        dplyr::mutate(prcg$bdata, dataset='pure')
    ) %>% dplyr::mutate(dataset = stringr::str_to_title(dataset))
    
    ggplot(hcounts, aes(x=dataset, y=hmrs, color=dataset, fill=dataset)) + 
        geom_violin(alpha=0.5) + 
        geom_jitter(width=0.2, shape=21, color='black') + 
        theme_bw() +
        labs(color="", fill="", x="", y="Hypomethylated regions from MethlySeekR\n(100bp windows)") + 
        theme_bw() + 
        scale_y_log10(n.breaks=10, labels=label_number(scale_cut = cut_si(''))) + 
        scale_color_manual(values=c('#BFBFBF', '#AAD2FF')) + 
        scale_fill_manual(values=c('#BFBFBF', '#AAD2FF')) + 
        theme(legend.position='none')
}

plot_cgi_counts <- function(scgi){
ggplot(scgi, aes(x=dataset, y=prop_osea, color=dataset)) + geom_boxplot() + geom_jitter(width=0.2) + 
    theme_minimal() + scale_y_continuous(limits=c(0,1)) + scale_color_manual(values=c("#BFBFBF", "#AAD2FF")) + 
    theme(legend.position='none') + labs(y="Proportion of HMRs in Open Seas", x="")
}

get_cdif_meth <- function(x){
    cdif = GRanges(seqnames='chr1', ranges=IRanges(start=156834317, end=160189186))
    cdol = findOverlaps(x$rHMRs, cdif)
    # Loci of rHMRs in the cd_if region
    cdix = x$rHMRs[queryHits(cdol)]$rHMRix
    cdhrc = x$rHMRc[rHMRix %in% cdix,]
    return(cdhrc)
}

plot_cdif <- function(dcdif, qcpass){
    if("meth" %in% names(dcdif)){
        dcdif$mean_meth=dcdif$meth
    }
    if("cg_ix" %in% names(dcdif)){
        dcdif$rHMRix = dcdif$cg_ix
    }
    log_info("Plotting p1")
    scdif = dcdif %>% dplyr::group_by(variable, dataset) %>% dplyr::summarise(mean_meth=mean(mean_meth, na.rm=T))
    p1_cdif = ggplot(scdif, aes(x=dataset, y=mean_meth, fill=dataset)) + geom_boxplot() + geom_jitter(alpha=0.5) # No major difference

    log_info("Plot P2")
    p2_cdif = ggplot(dcdif, aes(x=factor(rHMRix), y=mean_meth, color=dataset)) +  geom_smooth(se=FALSE, aes(group=variable), alpha=0.5) + facet_wrap(~dataset, scales='free_x')
    
    
    log_info("Plot P3")
    # Repeat the above except we give the Y axis as a specific sample order
    pwide = dcdif %>% dplyr::filter(dataset=='pure') %>% dplyr::select(rHMRix, variable, mean_meth) %>% 
        tidyr::pivot_wider(names_from=variable, values_from=mean_meth)
    pmat = dplyr::select(pwide, -rHMRix) %>% as.matrix() 
    # Cluster pmat by rows and then return the dendogram order of the column names
    porder = hclust(dist(t(pmat)))$order
    pnamesord = colnames(pmat)[porder]
    pnamesord = pnamesord[pnamesord %in% qcpass]
    dcord = dplyr::filter(dcdif, variable %in% pnamesord) %>% dplyr::mutate(vfact = factor(variable, levels=pnamesord))
    p3_cdif = ggplot(dcord, aes(x=factor(rHMRix), y=vfact)) + geom_tile(aes(fill=mean_meth)) + facet_wrap(~dataset, scales='free_x') + scale_fill_gradient2(low='blue', mid='white', high='red') + theme(axis.text.x=element_blank(), strip.background=element_blank()) + labs(x="CpGs", y="")

    if("immune" %in% names(dcdif)){
        log_info("Plot P4")
    icdif = dcdif %>% dplyr::group_by(variable, cls_immune, dataset, immune) %>% dplyr::summarise(mean_meth=mean(mean_meth, na.rm=T)) %>% dplyr::filter(!is.na(immune))
    p4_cdif = ggplot(icdif, aes(x=dataset, y=mean_meth, color=paste0(immune, cls_immune))) + geom_boxplot() + geom_jitter(alpha=0.5, position=position_jitterdodge(jitter.width=0.15)) + scale_color_manual(values=c("skyblue", "grey", "darkred", 'purple')) +
    theme_bw()
    }

    if("immune" %in% names(dcdif)){
        pall_cdif = patchwork::wrap_plots(p1_cdif, p2_cdif, p3_cdif, p4_cdif, ncol=1)
    }else{
    pall_cdif = patchwork::wrap_plots(p1_cdif, p2_cdif, p3_cdif, ncol=1)
    }
    return(pall_cdif)
}


temp_add_genes <- function(){

    eg = data.table::fread('data/source/cd1ifi6_regions.csv') %>%
        dplyr::rename(gene=hgnc_symbol, molecule=chromosome_name, start=start_position, end=end_position) %>%
        dplyr::mutate(orientation=1)
    eg = eg[40:80] %>% dplyr::filter(!stringr::str_starts(gene, "OR") & !gene=="")
    test_df = ggplot(dcord %>% dplyr::filter(dataset=='pure'), aes(x=factor(rHMRix), y=vfact)) + geom_tile(aes(fill=mean_meth)) + scale_fill_gradient2(low='blue', mid='white', high='red') + theme(axis.text.x=element_blank(), strip.background=element_blank()) + labs(x="CpGs", y="")
    gene_fig = ggplot(eg, aes(xmin = start, xmax = end, y = molecule, fill = gene, label=gene)) +
  geom_gene_arrow() +
    geom_gene_label(align = "left") + theme_bw() + theme(legend.position='none') + 
    scale_x_continuous(limits=c(min(eg$start), max(eg$end)))

    test_patch = patchwork::wrap_plots(test_df, gene_fig, ncol=1, heights=c(1, 0.1))
    ggsave('test.png', test_patch)

}

flatten_metho <- function(metho, dataset='bulk'){
      if(dataset=='bulk'){
    d_metho = list(
      pos=metho$pos, meth=metho$meth, cov=metho$cov
    )

     }else{
    d_metho = list(
      pos=metho$pos, meth=metho$meth_pure, cov=metho$cov_pure
    )

     }
    bmdata = cbind(d_metho$pos, d_metho$meth)
    bmdata$cg_ix = 1:nrow(bmdata)
    # using data.table, pivot to long format keeping chrom, start, end and cg_ix 
    bmdata = melt(bmdata, id.vars=c('chrom', 'start', 'end', 'cg_ix'))
    bmdata$dataset = dataset

    # Do as above but for d_metho$cov
    bcdata = cbind(d_metho$pos, d_metho$cov)
    bcdata$cg_ix = 1:nrow(bcdata)
    bcdata = melt(bcdata, id.vars=c('chrom', 'start', 'end', 'cg_ix'))
    bcdata$dataset = dataset

    # merge the two datasets
    setnames(bmdata, 'value', 'meth')
    setnames(bcdata, 'value', 'cov')
    bmdata = merge(bmdata, bcdata, by=c('chrom', 'start', 'end', 'cg_ix', 'variable', 'dataset'))
    return(bmdata)
}

get_dscores <- function(){
dscores = qs::qread('/nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/mcrpc_rna/data/mcrcp_danaher_scores.qs')

dclust = data.table(variable=names(dscores$clusters),cls_immune= dscores$clusters) %>% dplyr::mutate(
    immune=dplyr::case_when(
        cls_immune==1 ~ 'intermediate',
        cls_immune==2 ~ 'cold',
        cls_immune==3 ~ 'hot',
        cls_immune==4 ~ 'intermediate'
    )
)
return(dclust)
}
