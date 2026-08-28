# See scripts/main/dmr_lung_hilobulkpure.R

get_hilo_purity_tx <- function(s3, inv) {
    # Query patients with purity gaps (10)
    pgaps = s3 %>% dplyr::group_by(publication_id) %>%
        dplyr::summarise(rmin=min(rrbs_rho), rmax=max(rrbs_rho)) %>%
        dplyr::mutate(rdiff=rmax-rmin) %>%
        dplyr::arrange(-rdiff) %>% 
        dplyr::filter(abs(rdiff) > 0.25)
    pmins = s3 %>% dplyr::select(publication_id, region_min=region_id, rmin=rrbs_rho, pmin=rrbs_psi)
    pmaxs = s3 %>% dplyr::select(publication_id, region_max=region_id, rmax=rrbs_rho, pmax=rrbs_psi)
    pgaps = pgaps %>%
        left_join(pmins) %>%
        left_join(pmaxs)
    # Annotate LTX ID
    inv_ltx = dplyr::select(inv, publication_id, tracerx_id) %>% unique()
    pgaps = pgaps %>% left_join(inv_ltx)
    
    return(pgaps)
}

get_dmr_set <- function(x){ x %>% 
    dplyr::select(names(x)[!duplicated(names(x))]) %>% # Fix repeated entries if 'start' 
    dplyr::filter(!is.na(DMR)) %>% pull(cluster_id) %>% unique()}

count_dmrs <- function(x){ sum(!is.na(x$DMR))}

process_dmr_data <- function(metad, dmr_dir) {
    txidc = metad$tracerx_id %>% gsub('LTX0', 'LTX', .)
    shigh = metad$region_max
    slow = metad$region_min
    
    # Load DMR data for both samples
    hidir = path(dmr_dir, txidc, "Methylation", shigh, 're_segmented')
    lodir = path(dmr_dir, txidc, "Methylation", slow, 're_segmented')

    load(path(hidir, "bulk_tumour_promoter_DMR_calls.RData"))
    bhiprom = df_annotated_regions
    load(path(lodir, "bulk_tumour_promoter_DMR_calls.RData"))
    bloprom = df_annotated_regions
    load(path(hidir, "purified_tumour_promoter_DMR_calls.RData"))
    phiprom = df_annotated_regions
    load(path(lodir, "purified_tumour_promoter_DMR_calls.RData"))
    ploprom = df_annotated_regions

    # Create R object with DMR counts
    dmr_data = list(
        patient = txidc,
        metad = metad,
        dcount = list(
            bhiprom = count_dmrs(bhiprom),
            bloprom = count_dmrs(bloprom),
            phiprom = count_dmrs(phiprom),
            ploprom = count_dmrs(ploprom)
        ),
        setlist = list(
            bhiprom = get_dmr_set(bhiprom),
            bloprom = get_dmr_set(bloprom),
            phiprom = get_dmr_set(phiprom),
            ploprom = get_dmr_set(ploprom)
        ),
        setnames=c(
            "bhiprom"=metad_to_clean_name(metad, hi=T, pure=F),
            "bloprom"=metad_to_clean_name(metad, hi=F, pure=F),
            "phiprom"=metad_to_clean_name(metad, hi=T, pure=T),
            "ploprom"=metad_to_clean_name(metad, hi=F, pure=T),
            "patient"=metad$publication_id
        )
    )
    
    return(dmr_data)
}

metad_to_clean_name <- function(metad, hi=T, pure=T){
    # Get sample name for pretty printing in venn diagram
    reg = ifelse(hi, metad$region_max, metad$region_min)
    purity=round(ifelse(hi, metad$rmax, metad$rmin), 2)
    ploidy=round(ifelse(hi, metad$pmax, metad$pmin), 2)
    prefix=ifelse(hi, "Bulk", "Bulk (Lower Purity)")
    prefix=ifelse(pure, 'CAMDAC-purified', prefix)
    return((paste0(reg, ':', prefix, '\n', '(\u03C1=',purity,';\u03C8=',ploidy,')')))
}


plot_dmr_venn <- function(dmr_data) {
    # Setup data for plotting
    set_lists <- dmr_data$setlist[c("bhiprom", "bloprom", "ploprom")]
    set_names <- dmr_data$setnames[c("bhiprom", "bloprom", "ploprom")]
    names(set_lists) <- set_names
    pid <- dmr_data$setnames[['patient']]

    # Plot DMR venn diagram
    dmrvenn <- ggvenn(set_lists, text_size=4, set_name_size=5) + 
        labs(title=pid) + 
        scale_fill_manual(values=c('#865faa', '#a0a2a5', '#5aa3f4')) +
        theme(plot.title=element_text(hjust=0.5, vjust=-1)) + 
        scale_y_continuous(expand=c(0.1,0.1)) +
        scale_x_continuous(expand=c(0.1,0.1))
    
    return(dmrvenn)
}

get_er_stats <- function(dmr_dmata){
    bhi=dmr_data$setlist$bhiprom # in high bulk
    plo=dmr_data$setlist$ploprom # in camdac lo
    blo=dmr_data$setlist$bloprom # in bulk lo
    bhi_plo_ol = intersect(bhi, plo) # in high bulk and camdac lo 
    plo_enrich = setdiff(bhi_plo_ol, blo) # in high bulk and camdac lo, but not in bulk lo
    all_dmrs=unique(c(bhi, plo, blo))
    # as above, but not in bulk lo

    ol=intersect(bhi, plo)
    data.frame(
        patient=dmr_data$patient,
        bhi=length(bhi),
        plo=length(plo),
        blo=length(blo),
        plo_enrich=length(plo_enrich),
        perc_ploe_bhi=round(length(plo_enrich)/length(bhi)*100, 2),
        perc_ploe_all=round(length(plo_enrich)/length(all_dmrs)*100, 2)
    )
}
