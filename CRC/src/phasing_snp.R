parse_norm <- function(x){
    dnorm = fread(x)
    dcg = dnorm[width==2]
    dcg = dcg[!duplicated(dcg[, .(chrom, start, end)])]
    dcg = dcg[!is.na(m)]
    return(dcg)
}

# Boxplot as above, but add the N above each box
create_icr_boxplot_count <- function(dicd, fname) {
    bplt = ggplot(dicd, aes(x=id, y=m)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=3)) + geom_text(aes(x=id, label=..count..),size=2, y=1.0,stat='count', vjust=-1, inherit.aes=F) + scale_y_continuous(limits=c(0,1.5))
    ggsave(fs::path(figdir, fname), bplt, width=10, height=5)
}

# ICR
create_icr_median_hist <- function(dicd, fname) {
    dicmed = dplyr::group_by(dicd, id) %>% dplyr::summarise(med=median(m), cov=median(total_counts_m), ncg=n())
    dmeds = ggplot(dicmed, aes(x=med)) + geom_histogram(binwidth=0.05, color='skyblue', fill='skyblue') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() 
    dmeds_hc = ggplot(dicmed %>% dplyr::filter(cov>=5), aes(x=med)) + geom_histogram(binwidth=0.05, color='purple', fill='purple') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() 
    ggsave(fs::path(figdir, fname), wrap_plots(dmeds, dmeds_hc, ncol=1), width=5, height=10)
}

get_icr_med_window <- function(dicd){
    dicmed = dplyr::group_by(dicd, id) %>% dplyr::summarise(med=median(m), cov=median(total_counts_m), ncg=n())
    icr2 = dplyr::filter(dicmed, dplyr::between(med, 0.2, 0.8))
    dcln = merge(dicd, icr2, by='id')
    return(dcln)
}

helper_load_asm <- function(phase_dir, chrom="all", cln_output=T){
    # Load SNV phased data.
    # Returns a complete alignment dataset object.
    # This is required as the SNP (snv) POS is not in the ASM result, so we need the full alignment
    # To filter out SNVs at the CpG loci.
    # phase_dir: The directory containing the SNV phase results. See snv_mt_out
    # chrom: The chromosome to load. Default is 1. Set to "all" to load full dataset.
    
    shap = fs::dir_ls(phase_dir, glob="*hap_stats.csv.gz")
    smap = fs::dir_ls(phase_dir, glob="*phase_map.csv.gz")
    sinp = fs::dir_ls(phase_dir, glob="*snps.txt")
    stopifnot(all(fs::file_exists(c(shap, smap, sinp))))

    scna = fs::dir_ls(phase_dir, glob="*counts_cna.csv.gz")
    if(length(scna)==0){
        scna = fs::dir_ls(phase_dir, glob="*counts.csv.gz")
    }
    stopifnot(length(scna)>0)

    chrom_filter <- function(x, chr){
        if(chr=="all"){
            return(x)
        }else{
            return(x[x$chrom %in% chr])
        }
    }

    # Read data tables
    dt = fread_chrom(scna) %>% chrom_filter(chrom)
    dhap = fread_chrom(shap) %>% chrom_filter(chrom)
    dmap = fread_chrom(smap) %>% chrom_filter(chrom)
    dsip = fread_chrom(sinp) %>% chrom_filter(chrom) %>% dplyr::select(-BAF)
    
    # Merge data
    names(dsip)[2:ncol(dsip)] = paste0("hap_", names(dsip)[2:ncol(dsip)])
    names(dsip)[2] = "hap_POS" # fix for merge
    daln = merge(dhap, dmap, by=c("chrom", "hap_id"))
    daln = merge(daln, dsip, by=c("chrom", "hap_POS", "hap_ref", "hap_alt"), all.x=T)
    daln = merge(daln, dt, by=c("chrom", "start", "end"))


    if('cna_end' %in% names(daln)){
        extras = c("cna_start", "cna_end", "major_cn", "minor_cn", "purity", "ploidy")
    }else{
        extras = character()
    }

    if(cln_output){
        sel_names <- unique(c("chrom", "start", "end", "width", names(dsip),
            "hap_id", "hap_ref", "hap_alt", "hap_POS", "hap_BAF", "hap_reads", "hap_unexp",
            "alt_m", "alt_total_counts_m", "ref_m", "ref_total_counts_m", "qname", extras))
        dalnc = daln[, ..sel_names]
        setnames(dalnc, "alt_m", "altm")
        setnames(dalnc, "alt_total_counts_m", "atcm")
        setnames(dalnc, "ref_m", "refm")
        setnames(dalnc, "ref_total_counts_m", "rtcm")
        return(dalnc)
    } else {
        return(daln)
    }

}

fread_chrom <- function(x, ...) {
  # Read a file with a chrom column and ensure it is a factor
  # Used for counts and snp files where X and Y are missing
  x <- data.table::fread(x, ...)
  x$chrom <- as.character(x$chrom)
  if (stringr::str_starts(x$chrom[[1]], "chr")) {
    x$chrom <- gsub("chr", "", x$chrom)
  }
  return(x)
}

calculate_asm_statistics <- function(asm_cln) {
    asm_stat = asm_cln %>% 
        dplyr::group_by(id) %>% 
        dplyr::summarise(ncg=n(), n25=sum(abs(altm-refm)>=0.25), n50=sum(abs(altm-refm)>=0.5)) %>% 
        dplyr::mutate(p25=n25/ncg, p50=n50/ncg) %>% 
        dplyr::arrange(-ncg) 
    return(asm_stat)
}

create_icr_plot <- function(data) {
  ggplot(data, aes(x=id, y=altm-refm)) + 
    geom_boxplot() + 
    geom_jitter(alpha=0.2) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=3)) + 
    geom_text(aes(x=id, label=..count..), size=2, y=1.0, stat='count', vjust=-1, inherit.aes=F) + 
    scale_y_continuous(limits=c(-1,1)) + 
    theme_bw() 
}