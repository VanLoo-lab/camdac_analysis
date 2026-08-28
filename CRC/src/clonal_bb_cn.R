load_clonal_bb_47 <- function(bb_dir) {
  # Meta
  sample_id <- fs::path_file(bb_dir)

  # Load purity and ploidy
  pp_file <- fs::dir_ls(paste0(bb_dir, "/M-CallSubclones/"), glob = "*purity_ploidy.txt")
  if(length(pp_file) == 0){
    return(data.table())
  }
  pp <- data.table::fread(pp_file)
  fit_file <- fs::dir_ls(paste0(bb_dir, "/L-FitCopyNumber/"), glob = "*rho_and_psi.txt")
  if(length(fit_file) == 0){
    return(data.table())
  }
  ff <- suppressWarnings(data.table::fread(fit_file))
  pp$fit <- round(ff[is.best == T]$distance * 100, 3)

  # Load cna file
  cna_file <- fs::dir_ls(paste0(bb_dir, "/M-CallSubclones/"), glob = "*_subclones.txt")
  cna_file <- cna_file[!grepl("chrX", cna_file)]
  bb_cna <- data.table::fread(cna_file)
  bb_cna_fields <- c(
    "chr", "startpos", "endpos", "ntot", "nMaj1_A", "nMin1_A",
    "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A"
  )
  bb_cna <- bb_cna[, ..bb_cna_fields]
  bb_cna <- bb_cna[!is.na(nMaj1_A)]

  # Ensure Maj/Min2 are integer. Required when NA in field for compleetely clonal sample is read
  bb_cna$nMaj2_A <- as.integer(bb_cna$nMaj2_A)
  bb_cna$nMin2_A <- as.integer(bb_cna$nMin2_A)

  # Select major copy number from joint solutions
  bb_cna[, selector := (bb_cna$frac1_A > bb_cna$frac2_A) | is.na(bb_cna$frac2_A)]
  bb_cna[, nA := data.table::fifelse(
    selector, nMaj1_A, nMaj2_A
  )]
  bb_cna[, nB := data.table::fifelse(
    selector, nMin1_A, nMin2_A
  )]
  bb_cna[, CN := nA + nB]

  # Finalise data for export
  cna_clean <- bb_cna[, .(
    chrom = factor(chr, levels = c(1:22, "X", "Y")), start = startpos, end = endpos,
    nA, nB, CN
  )]
  setkeyv(cna_clean, cols = c("chrom", "start", "end"))

  cna_clean$sample_id <- sample_id
  cna_clean$purity <- pp$purity
  cna_clean$ploidy <- pp$ploidy
  cna_clean$fit <- pp$fit

  # For CAMDAC paper: add fraction of genome LOH
  bp_loh = cna_clean[ nB==0, .(bp=sum(end-start)) ]$bp/1e4 # Convert to Mb to avoid int. overflow
  bp_het = cna_clean[ nB>0, .(bp=sum(end-start)) ]$bp/1e4 # Convert to Mb to avoid int. overflow
  cna_clean$frac_loh = bp_loh / (bp_loh + bp_het)

  return(cna_clean)
}


aggregate_battenberg_cna_segments_47 <- function(bb_dirs, outfile){
    # Assumes all folders in resdirs contain a single sample's battenberg output
    # Returns a table with CNA calls for all samples.
    bb_cna = lapply(bb_dirs, load_clonal_bb_47)
    bb_cna = data.table::rbindlist(bb_cna)
    data.table::fwrite(bb_cna, outfile)
    print(sprintf("Aggregated CNA segments saved to: %s", outfile))
}