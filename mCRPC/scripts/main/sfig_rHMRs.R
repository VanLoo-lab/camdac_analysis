# Generate figures from rHMRs replicating Zhao et al 2020

library(data.table)
library(dplyr)
library(ggplot2)
library(logger)
library(ComplexHeatmap)
options(bitmapType='cairo')

# Select hmr files
bhmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/bulk_hmr.qs'
phmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/pure_hmr.qs'
thmr = 'data/out/hmrs/test_bulk_res.qs'
meta = fread('data/out/hmrs/meta_mcrpc.tsv')

# Run HMR analysis pipeline

hmr_fig_pipe <- function(hmr_file, meta, fig_outfile, outfile){

    if(fs::file_exists(outfile)){
        log_info("Skipping {outfile}")
        return(NA)
    }

    log_info("START:::{outfile}")
    fs::dir_create(fs::path_dir(fig_outfile))
    fs::dir_create(fs::path_dir(outfile))
    
    # Data for 1 HMR cohort (bulk or pure)
    log_info("Reading HMR data")
    data = qs::qread(hmr_file)

    ## Transformation for later. Format count of cpgs iwthin hmrs
    ccount = data$rHMRc[, .(ncg_meth_inhmrs=sum(ncg_meth, na.rm=T)), by='variable']
    names(ccount) = c("sample_id", "ncg_meth_inhmrs")

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

    # Get the top N% of rHMRs with data in at least half of the samples
    half_samps = round(ncol(rHMRwide[, -1])/2, 0)
    rvd = data$rHMRv[nsamp >= half_samps,]
    rvd = rvd[order(-rvd$rvar),]
    top10pc = round(nrow(rvd)*.10,0)
    topix = rvd[1:top10pc]$rHMRix # Expectd to be sorted by rvar
    rtops = rHMRwide[rHMRix %in% topix,]
    rtops$rHMRix = NULL
    rtops_mat = as.matrix(rtops)
    rtops_mat[is.na(rtops_mat)] = 0.5

    # Setup annotations based on the metadata table
    rtmeta = meta[match(names(rtops), meta$sample_id),]
    rtmeta_s = dplyr::select(rtmeta, -sample_id)


    # Create complex heatmap from rtops using data$hmrpat as annotation
    log_info("Building heatmap")
    hm = Heatmap(
        rtops_mat,
        name = "Methylation",
        column_title = "Samples",
        row_title = "rHMRs",
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = HeatmapAnnotation(df = rtmeta_s),
        heatmap_legend_param = list(title = "M", at = c(0, 0.5, 1))
    )

    # Next, generate the HMR boxplots
    log_info("Calculating boxplots")
    bdata = merge(data$hmatc, meta, by='sample_id', all.x=T)
    bdata = merge(bdata, data$hmr_ps, by='sample_id', all.x=T)
    bdata = merge(bdata, ccount, by='sample_id', all.x=T)
    bdata = dplyr::mutate(bdata,
        CSTATUS=ifelse(cimp==T, "CIMP", "Non-CIMP"))
    bdata$CSTATUS = factor(bdata$CSTATUS, levels=c("Non-CIMP", "CIMP"))

    h1 = ggplot(bdata, aes(x=CSTATUS, y=hmrs)) +
        geom_boxplot() +
        geom_jitter(aes(color=CSTATUS)) +
        theme_classic() +
        scale_color_manual(values=c("black", "#72bcd4")) +
        labs(y="HMRs per sample (thousands)", x="") +
        theme(legend.position="none")

    h2 = ggplot(bdata, aes(x=CSTATUS, y=ncg_meth_inhmrs)) +
        geom_boxplot() +
        geom_jitter(aes(color=CSTATUS)) +
        theme_classic() +
        scale_color_manual(values=c("black", "#72bcd4")) +
        labs(y="CpGs methylated within rHMRs", x="") +
        theme(legend.position="none")
        

        ## TODO: Shown as a percentage in the paper?

    h3 = ggplot(bdata, aes(x=CSTATUS, y=i_s_s/1000)) +
        geom_boxplot() +
        geom_jitter(aes(color=CSTATUS)) +
        theme_classic() +
        scale_color_manual(values=c("black", "#72bcd4")) +
        labs(y="rHMRs in CpG islands, shores and shelves (thousands)", x="") +
        theme(legend.position="none")

        ## TODO: Article has error. rHMRs not HMRs

    h4 = ggplot(bdata, aes(x=CSTATUS, y=opensea/1000)) +
        geom_boxplot() +
        geom_jitter(aes(color=CSTATUS)) +
        theme_classic() +
        scale_color_manual(values=c("black", "#72bcd4")) +
        labs(y="rHMRs in open seas (thousands)", x="") +
        theme(legend.position="none")

        ## TODO: Article has error. rHMRs not HMRs

    # Combine and save figure
    log_info("Saving figure")
    hmrplots = patchwork::wrap_plots(h1, h2, h3, h4, ncol=2)
    heatplot = grob = grid.grabExpr(draw(hm)) 
    ggsave(fig_outfile, patchwork::wrap_plots(heatplot, hmrplots, nrow=1), width=25, height=15, dpi=300)

    # Save object to the output results file
    log_info("Saving results")
    qs::qsave(list(rtops, bdata, hmrplots, heatplot), outfile)

    log_info("Complete")
}

# Run the pipeline
hmr_fig_pipe(thmr, meta, 'data/out/hmrs/test_hmr_fig.png', 'data/out/hmrs/res_hmr_fig.qs')
hmr_fig_pipe(bhmr, meta, 'data/out/hmrs/bulk_hmr_fig.png', 'data/out/hmrs/res_bulk_hmr_fig.qs')
hmr_fig_pipe(phmr, meta, 'data/out/hmrs/pure_hmr_fig.png', 'data/out/hmrs/res_pure_hmr_fig.qs')
