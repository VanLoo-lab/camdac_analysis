library(GenomicRanges)
library(data.table)
library(janitor)
library(dplyr)
library(foreach)
library(doParallel)
library(logger)

load_wang_data <- function(){
    wang <- data.table::fread('data/source/luad_lusc_wang_2020.csv')
    gregions <- data.table::fread('data/source/240416_luad_lusc_atac_rna_regions.tsv') %>% janitor::clean_names() %>%
        dplyr::mutate(gene=gene_name)
    data <- merge(wang, gregions, by='gene')
    data$rid = seq(1, nrow(data))
    return(data)
}

load_wang_enhancers <- function(){
    data.table::fread('/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/at2_basal_enhancers.csv')
}

ranges_and_proms <- function(data, enh){
  # Create a GRanges object for gene regions
  genes = GRanges(seqnames = data$chromosome_scaffold_name,
                   ranges = IRanges(start = data$gene_start_bp, end = data$gene_end_bp),
                   strand = data$strand,
                   rid = data$rid,
                   gene = data$gene_name,
                   cell = data$cell,
                   dset = 'genes'
                   )

  # Add promoter (Using 1.5kb as per rest of CAMDAC paper)
  proms = promoters(genes, upstream=1500, downstream=0)
  proms$dset = 'promoters'

  # Renh
  names(enh)[1] = 'gene'
  enh = merge(enh, unique(data[, .(gene, cell)]), all.x=T, by='gene')
  renh = GRanges(seqnames = enh$chrom,
                  ranges = IRanges(start = enh$start, end = enh$end),
                  strand = '*',
                  rid = seq(1, nrow(enh)),
                  gene = enh$gene,
                  cell = enh$cell,
                  dset = 'enhancers'
  )

  # Combine
  areg = c(genes, proms)
  areg = c(areg, renh)

  # Add enhancer
  return(areg)
}


load_sample_anno <- function(){
 inv = data.table::fread('/camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/post_inventory_qc_ascat_22_05_19.tsv') %>%
     dplyr::mutate(pid=gsub("LTX0", "LTX", tracerx_id)) %>%
     dplyr::select(cruk=TRACERxID, pid, region_id)
 clin = data.table::fread('/camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/supp_clinical_luad_lusc_candisc.csv') %>%
     dplyr::mutate(subtype=ifelse(Diagnosis=="Invasive adenocarcinoma", "LUAD", "LUSC")) %>%
     dplyr::select(cruk=Publication_ID, subtype)
 keep = merge(inv, clin, by='cruk')
 return(keep)
}

parse_sdata <- function(ff, sanno){
    intype = data.table(file=ff) %>%
        dplyr::mutate(pid = stringr::str_extract(file, 'LTX\\d+')) %>%
        dplyr::mutate(region_id = file %>% stringr::str_extract('Methylation/[A-Za-z]+\\d+') %>% gsub('Methylation/', '', .)) 
    intype = merge(intype, sanno, by=c("pid", "region_id"))
    intype$file=NULL
    return(intype)
}


process_meth_dmp <- function(x, areg){
    # Get GRANGES from x
    gr = GRanges(seqnames = x$chromosome, ranges = IRanges(start = x$start, end = x$end), strand = '*')
    # For each element of areg, get all the indexes of x that overlap
    idx = findOverlaps(areg, gr)
    oldata = data.table(cbind(
        areg[queryHits(idx)],
        x[subjectHits(idx),]
    ))
    # Summarise the dataset 
    olsum = oldata[, .(
        mt_av = mean(m_t_corr),
        mn_av = mean(m_n),
        mb_av = mean(m_b),
        mdiff_av = mean(m_diff),
        prob_hyper_av = mean(prob_hyper, na.rm=T),
        nCpGs = .N,
        DMP_t_frac = sum(!is.na(DMP_t_corr))/.N,
        DMP_t_hy = sum(DMP_t_corr=="hyper" & !is.na(DMP_t_corr)),
        DMP_t_ho = sum(DMP_t_corr=="hypo" & !is.na(DMP_t_corr)),
        DMP_t_non = sum(is.na(DMP_t_corr)),
        DMP_b_frac = sum(!is.na(DMP_b_corr))/.N,
        DMP_b_hy = sum(DMP_b_corr=="hyper" & !is.na(DMP_b_corr)),
        DMP_b_ho = sum(DMP_b_corr=="hypo" & !is.na(DMP_b_corr)),
        DMP_b_non = sum(is.na(DMP_b_corr)),
        cov_b_av = mean(cov_b),
        cov_t_av = mean(cov_t)
    ), by=c('rid', 'dset', 'gene', 'cell')]
    return(olsum)
}