
calc_n_meth_pass <- function(meth_file, dir){
    outfile = path(dir, paste0(basename(meth_file), ".n_meth_pass.tsv"))
    if(file.exists(outfile)){
        return(outfile)
    }
    # Load data
    dt = data.table::fread(meth_file, select=c("chrom", "start", "end", "m", "cov"))[!is.na(m)]
    # Calculate number of CpGs with >5x coverage
    dt_result = dt %>% 
        summarise(
            n1= sum(cov >= 1),
            n5= sum(cov >= 5),
            n10= sum(cov >= 10),
            n30= sum(cov >= 30),
            n100= sum(cov >= 100),
            ncgs=n()
        ) %>% 
        mutate(
            n5_perc = sum(n5*100/ncgs),
            n10_perc = sum(n10*100/ncgs),
            n30_perc = sum(n30*100/ncgs),
            n100_perc = sum(n100*100/ncgs)
        ) %>% 
        mutate(fname = basename(meth_file) %>% gsub(".m.csv.gz", "", .))
    
    fwrite(dt_result, outfile, sep="\t")
    return(outfile)
}

calc_cov_boxplot <- function(ac_file, dir){ 

  outfile = path(dir, paste0(basename(ac_file), ".cov_boxplot.tsv"))
  if(file.exists(outfile)){
    return((outfile))
  }
  # Load data
  dt = data.table::fread(ac_file, select=c("total_counts", "total_depth", "m", "POS"))[!is.na(m)]
  
  # Calculate tumour coverage
  dt_result = dt %>% 
    summarise(total_counts= list( boxplot.stats(total_counts)$stats),
              total_depth=list(boxplot.stats(total_depth)$stats),
              stats=list( c('lower_whisker','lower_hinge','median','upper_hinge','upper_whisker') ),
              n_meth_n1 = n()) %>% 
    unnest(c(total_counts, total_depth, stats)) %>% 
    pivot_longer(c(total_counts, total_depth), names_to="count_type") %>% spread(stats, value)

  dt_result = cbind(
    data.frame(fname=basename(ac_file) %>% gsub(".SNPs.CpGs.all.sorted.csv.gz", "", .)),
    dt_result
  )

  fwrite(dt_result, outfile, sep="\t")
  return(outfile)
}


# Table script for mcrpc

calc_nrpcc <- function(cov, purity, ploidy) {
  if(any(sapply(c(cov, purity, ploidy), is.na))){
    return(NA)
  }
  purity * cov /
    ((purity*ploidy) + ((1-purity) * 2))
}

v_calc_nrpcc <- Vectorize(calc_nrpcc)


su2c_id <- function(x) str_extract(x, "DTB-\\d+")
read_pp_data <- function(x) {
    fread(x) %>% mutate(patient=su2c_id(x))
}
read_snv_data <- function(x) {
    fread(x) %>% mutate(super=ccf>=1.5, patient=su2c_id(x)) %>% 
        dplyr::group_by(patient) %>% dplyr::summarise(
            count = n(),
            nsuper = sum(super),
            frac_super = nsuper/count,
        ) %>% as.data.table
}
read_tissue <- function(x){
    x %>% fread %>% filter(ANALYTE_TYPE=="DNA" & HISTOLOGICAL_TYPE!="white blood cells") %>%
    mutate(patient = su2c_id(SAMPLE_ID), IS_TUMOR=ifelse(IS_TUMOR=="yes", T, F)) %>%
    dplyr::select(patient, is_tumor=IS_TUMOR, body_site=BODY_SITE) %>%
    dplyr::mutate(body_site = ifelse(body_site=="LN", "lymph", body_site))
}

