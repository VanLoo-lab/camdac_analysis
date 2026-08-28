# Load dependencies
library(data.table)

# Load data
cgi_ucsc_hg38 = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz'
dt = data.table::fread(cgi_ucsc_hg38)

# Reformat data
dt = dt[, .(name=V1, chrom=V2, start=V3, end=V4, length=V4-V3+1)]
dt = dt[chrom %in% c(paste0("chr", 1:22), "chrX", "chrY")]
dt$chrom = gsub("chr", "", dt$chrom)
dt$cgix = 1:nrow(dt) # Index for each CGI
dt$dtype = "islands" # Label as CGIs

# Generate a new data table annotation for shores, shelves and open seas, as a function of island annotation
dt_sso = rbind(
    dt[, .(chrom, start=start-2000, end=start-1, name, length=2000, cgix, dtype="shore")], # Left shore
    dt[, .(chrom, start=end+1, end=end+2000, name, length=2000,  cgix, dtype="shore")], # Right shore
    dt[, .(chrom, start=start-4000, end=start-2001, name, length=2000, cgix, dtype="shelf")], # Left shelf
    dt[, .(chrom, start=end+2001, end=end+4000, name, length=2000,cgix,  dtype="shelf")], # Right shelf
    dt[, .(chrom, start=start-6000, end=start-4001, name, length=2000,cgix,  dtype="opensea")],
    dt[, .(chrom, start=end+4001, end=end+6000, name, length=2000, cgix, dtype="opensea")]
)

# Combine with islands
dt = rbind(dt, dt_sso)

# Arrange by genomic coordinates
setkey(dt, chrom, start, end)
dt = dt[order(chrom, start, end)]

# Write to output file
data.table::fwrite(dt, "data/meta/hg38_cgi_clean.csv.gz")
