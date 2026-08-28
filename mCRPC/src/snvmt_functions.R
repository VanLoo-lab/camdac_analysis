# Helper functions to convert CAMDAC ASM pipeline to SNV phasing
snv_mt_pipeline <- function(sample, config, ccf_cutoff=0.7) {
    # Create the results directory
    sout <- fs::dir_create(snv_mt_out(sample, config))

    # Check for expected output. If all exist, skip.
    outfiles <- fs::path(sout, paste0(
        paste0(sample$patient_id, ".", sample$id, "."),
        c(
            "snv_counts.csv.gz", "snv_hap_stats.csv.gz", "snv_phase_map.csv.gz",
            "snv_snps.txt", ".snv_counts_cna.csv.gz"
        )
    ))
    if (all(fs::file_exists(outfiles))) {
        loginfo("All SNV phase files exist. Skipping.")
        return(path(sout, "snv_counts.csv.gz"))
    }

    loginfo("Running SNV phase pipeline")
    
    apply_ccf_cutoff(sample, config, ccf_cutoff)

    # Run allele counting for SNV mt
    # 1. Adapt write_asm_counts_output to a temporary object that we can load later
    testthat::local_mocked_bindings(
        write_asm_counts_output = hack_write_asm_counts_output,
        .package = "CAMDAC"
    )
    csnv_temp <- cmain_asm_allele_counts(sample, config)
    # csnv_temp = path(sout, "temp_obj.qs")

    # Write results object to appropriate files
    loginfo("Writing SNV phase files")
    csnv <- qs::qread(csnv_temp)
    fwrite(csnv$asm_cg, outfiles[1])
    fwrite(csnv$hap_stats, outfiles[2])
    fwrite(csnv$map, outfiles[3])

    # Move the asm_snps file to the tsnv snps file
    asm_snps <- get_fpath(sample, config, "asm_snps")
    fs::file_move(asm_snps, outfiles[4])

    # Remove the temporary object
    fs::file_delete(csnv_temp)

    # Add CNA calls to the snv counts file
    loginfo("Annotating SNV counts file with CNA calls")
    anno_out <- path(sout, paste0(sample$patient_id, ".", sample$id, ".snv_counts_cna.csv.gz"))
    csnv_cna_file <- annotate_snv_mt_cna(outfiles[1], sample, config, anno_out)

    # Return
    loginfo("SNV phase pipeline complete")
    return(csnv_cna_file)
}

hack_write_asm_counts_output <- function(result, sample, config) {
    # sout code taken from outside function as not yet in get_fpath
    loginfo("HACK: Writing SNV counts temp output")
    sout <- fs::path(config$outdir, sample$patient_id, "SNV_phase", sample$id)
    temp_obj <- fs::path(sout, "temp_obj.qs")
    qs::qsave(result, temp_obj)
    return(temp_obj)
}

annotate_snv_mt_cna <- function(snv_counts_file, sample, config, anno_out) {
    if (fs::file_exists(anno_out)) {
        loginfo("Annotated SNV counts file exists. Skipping.")
        return(anno_out)
    }

    loginfo("Annotating SNV counts file with CNA calls")
    # Read counts
    scounts <- fread_chrom(snv_counts_file)

    # Read cna
    cna <- get_fpath(sample, config, "cna")
    cna <- fread(cna)

    # Annotate CNA
    setkey(scounts, chrom, start, end)
    setkey(cna, chrom, start, end)
    scover <- foverlaps(scounts, cna, type = "any", mult = "first", nomatch = NA)
    scover[, `:=`(
        cna_start = start,
        cna_end = end,
        start = i.start,
        end = i.end,
        i.start = NULL, # drop
        i.end = NULL # drop
    )]

    # Write to file
    fwrite(scover, anno_out)

    return(anno_out)
}

helper_load_snv_phased <- function(phase_dir, chrom="all", cln_output=T){
    # Load SNV phased data.
    # Returns a complete alignment dataset object.
    # This is required as the SNP (snv) POS is not in the ASM result, so we need the full alignment
    # To filter out SNVs at the CpG loci.
    # phase_dir: The directory containing the SNV phase results. See snv_mt_out
    # chrom: The chromosome to load. Default is 1. Set to "all" to load full dataset.
    
    scna = fs::dir_ls(phase_dir, glob="*.snv_counts_cna.csv.gz")
    shap = fs::dir_ls(phase_dir, glob="*.snv_hap_stats.csv.gz")
    smap = fs::dir_ls(phase_dir, glob="*.snv_phase_map.csv.gz")
    sinp = fs::dir_ls(phase_dir, glob="*.snv_snps.txt")
    stopifnot(all(fs::file_exists(c(scna, shap, smap, sinp))))

    dt = fread_chrom(scna)
    if(chrom!="all"){
        csel = chrom
        dt = dt[chrom %in% csel]
    }

    dhap = fread_chrom(shap)
    dmap = fread_chrom(smap)
    dsip = fread_chrom(sinp)
    # Set SNV names to contain haplotype before merge
    names(dsip)[2:ncol(dsip)] = paste0("hap_", names(dsip)[2:ncol(dsip)])
    names(dsip)[2] = "hap_POS" # fix for merge
    daln = merge(dhap, dmap, by=c("chrom", "hap_id"))
    daln = merge(daln, dt, by=c("chrom", "start", "end"))
    daln = merge(daln, dsip, by=c("chrom", "hap_POS", "hap_ref", "hap_alt"), all.x=T)

    if(cln_output){
        sel_names <- unique(c("chrom", "start", "end", "width", names(dsip),
            "hap_id", "hap_ref", "hap_alt", "hap_POS", "hap_BAF", "hap_reads", "hap_unexp",
            "alt_m", "alt_total_counts_m", "ref_m", "ref_total_counts_m",
            "cna_start", "cna_end", "major_cn", "minor_cn", "purity", "ploidy", "qname"))
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

snv_mt_out <- function(sample, config){
    fs::path(config$outdir, sample$patient_id, "SNV_phase", sample$id)
}

apply_ccf_cutoff <- function(sample, config, ccf_cutoff){
    # Cache input asm_snp file, limit to ccf cutoff for pipeline
    sout <- snv_mt_out(sample, config)
    asm_snps <- get_fpath(sample, config, "asm_snps")
    asm_snps_cached <- fs::path(sout, paste0(sample$patient_id, ".", sample$id, ".snv_snps_input.txt"))
    fs::file_copy(asm_snps, asm_snps_cached, overwrite = TRUE)
    snv_data <- fread(asm_snps)
    if(ncol(snv_data) > 4){
        log_info("Filtering ASM SNPs to CCF >= {ccf_cutoff}")
        # Column 5 must be CCFs if present.
        ccf_bool = snv_data[[5]] >= ccf_cutoff & snv_data[[5]] <= 1.1 & !is.na(snv_data[[5]])
        snv_filtered = snv_data[ccf_bool]
        log_info("SNV count before filtering: {nrow(snv_data)}")
        log_info("SNV count after filtering: {nrow(snv_filtered)}")
        fwrite(snv_filtered, asm_snps)
    } else {
        log_info("No CCF column found. Skipping filtering.")
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
