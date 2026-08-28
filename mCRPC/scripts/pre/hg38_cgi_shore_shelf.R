# Make CGI annotation
library(data.table)
cgi_ucsc_hg38 = 'data/meta/cpgIslandExt.hg38.bed'

# Read data
dt = fread(cgi_ucsc_hg38)[, .(
    chrom=V2,
    start=V3,
    end=V4,
    name=V1,
    length=V4-V3+1
)]
dt = dt[chrom %in% c(paste0("chr", 1:22), "chrX", "chrY")]
dt$dtype = "islands"
dt$ix = 1:nrow(dt)

# Create df for shores, shelves and open seas
dt_sso = rbind(
    dt[, .(chrom, start=start-2000, end=start-1, name, length=2000, ix, dtype="shore")], # Left shore
    dt[, .(chrom, start=end+1, end=end+2000, name, length=2000,  ix, dtype="shore")], # Right shore
    dt[, .(chrom, start=start-4000, end=start-2001, name, length=2000, ix, dtype="shelf")], # Left shelf
    dt[, .(chrom, start=end+2001, end=end+4000, name, length=2000,ix,  dtype="shelf")] # Right shelf
    #dt[, .(chrom, start=start-6000, end=start-4001, name, length=2000,ix,  dtype="opensea")],
    #dt[, .(chrom, start=end+4001, end=end+6000, name, length=2000, ix, dtype="opensea")]
)
dt = rbind(dt, dt_sso)
dt$all_ix = 1:nrow(dt)

setkey(dt, chrom, start, end)
dt = dt[order(chrom, start, end)]

dt$chrom = gsub("chr", "", dt$chrom)
dt$cgi_width = dt$end - dt$start + 1
dt$cgi_annot = dt$dtype
dt$cgi_ix = dt$ix
dt = dt[, .(
    chrom, start, end, cgi_width, cgi_annot, cgi_ix, all_ix
)]

# Add CAMDAC WGBS CpG loci, taking any non-overlapping loci to be open seas
loci_files = fs::dir_ls('/nemo/project/proj-vanloo/secure/working/nmensah/pipeline_files/wgbs/hg38/loci_files')

parse_loci_file <- function(lf){
    lo = load(lf)
    xv = loci_subset[width(loci_subset)==2]
    seqlevelsStyle(xv) = "NCBI"
    xv = data.table(data.frame(xv))
    xv = xv[, .(chrom=seqnames, start, end)]
    # Warn: CpGs in loci file incomplete
    setkey(xv, chrom, start, end)
    setkey(dt, chrom, start, end)
    # Find overlapping ranges
    fol = data.table::foverlaps(xv, dt, type="any", mult='first', which=TRUE)
    xbind = cbind(xv, dt[fol,])
    xbind = xbind[, .(chrom, start, end, cgi_annot)]
    xbind$cgi_annot = ifelse(is.na(xbind$cgi_annot), "opensea", xbind$cgi_annot)
    return(xbind)
}

library(parallel)
dparsed = mclapply(loci_files, parse_loci_file, mc.cores=10)
dparsed = data.table::rbindlist(dparsed)
data.table::fwrite(dparsed, "data/meta/hg38_cgi_open_cln.csv.gz")

