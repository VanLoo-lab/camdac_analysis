library(data.table)
library(ggplot2)
library(matrixStats)
indir <- fs::dir_create("data/out/panel")
valfiles <- fs::dir_ls(fs::path(indir, "validation"))

outdir <- fs::dir_create("data/out/panel/val_rec_cgi")

# Load meta for cg locations
hg19_cgi <- fst::read_fst("/nemo/project/proj-vanloo/secure/working/nmensah/pipeline_files/rrbs/pipeline_files/hg19_annotations/hg19_all_CpGs_annotations.fst")
hg19_cg <- data.table(hg19_cgi)[, .(chrom, start, end, CGI)]
# hg38 cgi created manually from sea data as the CGI field of the CAMDAC annotation has no shelf label.
hg38_cg <- data.table::fread("/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/source/hg38_cgi_ucsc.csv.gz")

# Setup function to calculate result from validation file
calculate_rec_cgi <- function(x, hgsel, fname) {
    # For each patient CpG, get a vector of the positions that were DMPs, setting 1 for hyper, -1 for hypo and 0 for not
    xd <- dplyr::select(x, contains("DMP_b"))
    xd[is.na(xd) | xd == "N.S."] <- 0
    xd[xd == "hyper"] <- 1
    xd[xd == "hypo"] <- -1

    # We will now classify each CpG into the following classes:
    #   rec_hyper, rec_hypo, rec_both, non_rec
    nsamples <- ncol(xd)
    r50 <- floor(nsamples / 2)
    resdt <- data.table(
        nsamples = nsamples,
        rec_hyper = rowSums(xd == 1) >= 3,
        rec_hypo = rowSums(xd == -1) >= 3,
        rec_both = rowSums(xd == 1 | xd == -1) >= 3,
        rec50_hyper = rowSums(xd == 1) >= r50,
        rec50_hypo = rowSums(xd == -1) >= r50,
        rec50_both = rowSums(xd == 1 | xd == -1) >= r50
    )
    resdt[, `:=`(
        rec_none = !rec_hyper & !rec_hypo & !rec_both,
        rec50_none = !rec50_hyper & !rec50_hypo & !rec50_both
    )]
    resdt$fname <- fname
    resdt <- cbind(x[, .(chrom, start, end)], resdt)

    # Finally, annotate cpg islands
    # This is delicate as lung needs the RRBS hg19 annotations while WGBS need hg38
    setkey(resdt, chrom, start, end)
    setkey(hgsel, chrom, start, end)
    res <- merge(resdt, hgsel, all.x = T, by = c("chrom", "start", "end"))

    res$CGI <- ifelse(is.na(res$CGI), "opensea", res$CGI)
    res
}

library(foreach)
library(doParallel)
registerDoParallel(5)
foreach(vfile=valfiles) %dopar% {
    print(vfile)
    vald <- qs::qread(vfile)$pos
    fname <- fs::path_ext_remove(fs::path_file(vfile))
    outfile = fs::path(outdir, paste0(fname, ".csv.gz"))
    if(fs::file_exists(outfile)) {
        return(NULL)
    }
    if (!grepl("normal", fname)) {
        return(NULL)
    }
    # Get the hg annotation
    if(grepl("lung", fname)){
        hgsel = hg19_cg
    }else{ hgsel= hg38_cg }
    
    res <- calculate_rec_cgi(vald, hgsel, fname)
    data.table::fwrite(res, outfile)
}
