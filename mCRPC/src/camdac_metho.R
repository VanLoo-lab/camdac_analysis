# Load CAMDAC DNAme files as methylation objects (list)
library(data.table)
library(dplyr)
library(logger)
library(foreach)
library(doParallel)
library(itertools)
library(stringr)
library(GenomicRanges)

load_camdac_metho <- function(manifest, chunksize=10, ncores=3){
    # TODO: TESTING
    manifest = manifest
    
    # manifest: data.frame pointing to files with beta values for loading CAMDAC methylation data
    i_samp = seq_len(nrow(manifest))
    i_samp_chunk = ichunk(i_samp, chunksize, mode='integer')

    # sample names must be unique
    uniq_names = length(unique(manifest$sample_id)) == nrow(manifest)

    # files must exist
    all_exist = all(sapply(manifest$ac_file, fs::file_exists))

    # Return empty data table if any of the above conditions are not met
    if (!uniq_names | !all_exist) {
        log_warn("Manifest does not meet requirements. Returning empty data.table")
        return(data.table())
    }

    # Create a temporary directory to hold objects
    tdir = tempdir()
    
    ## FOR EACH SAMPLE CHUNK
    doParallel::registerDoParallel(cores=ncores)
    metho_l = foreach(i_list = i_samp_chunk, .combine='c') %dopar% {
            reslist = lapply(
                i_list,
                function(x){
                    ac_file = manifest[x,]$ac_file
                    sample_id = manifest[x,]$sample_id
                    log_info("Processing {sample_id}")
                    metho = dispatch_read_metho(sample_id, ac_file)
                    return(metho)
                }
            )
            if(length(reslist)==1){
                metho = reslist[[1]]
            }else{
                metho = Reduce(merge_metho, reslist)
            }
            # create tmpfile to save objetc
            tmpfile = tempfile(tmpdir=tdir)
            qs::qsave(metho, tmpfile)
            return(tmpfile)
        }
    log_info("Combining {length(metho_l)} metho objects...")
    metho = qs::qread(metho_l[[1]])
    remainder = metho_l[-1]
    for(i in seq_along(remainder)){
        log_info("Merging {i+1}...")
        obj = qs::qread(remainder[[i]])
        metho = merge_metho(metho, obj)
    }
    log_info("Complete!")
    
    doParallel::stopImplicitCluster()

    return(metho)

}

dispatch_read_metho <- function(sid, x){
    if (str_ends(x, "CpGs.all.sorted.csv.gz")){
        return(read_normal_as_metho(sid, x))
    }
    if (str_ends(x, ".bg.gz")){
        return(read_bed_metho(sid, x))
    }
    if (str_ends(x, "pure.csv.gz")){
        return(read_pure_metho(sid, x))
    }
    if (str_ends(x, "m.RData")){
        return(read_rrbs_metho(sid, x))
    }

    log_warn("Could not determine file type for {x}")
    return(data.table())
}

read_rrbs_metho <- function(sid, x){
    load(x)
    d = gr_normal_m
    seqlevelsStyle(d) = "NCBI"
    d = d[!is.na(d$m_n)]
    metho = metho_make(as.character(seqnames(d)), start(d), end(d), d$m_n, d$M_n + d$UM_n, sid)
    metho$issnp = data.table(x=d$SNP_n!=0)
    names(metho$issnp) = sid
    return(metho)
}

read_bed_metho <- function(sid, x){
    bed = fread(x, verbose=FALSE, showProgress = FALSE)
    names(bed) = c("chrom", "start", "end", "m")
    bed$chrom = gsub("chr", "", bed$chrom)
    metho = metho_make(bed$chrom, bed$start, bed$end, bed$m, 50, sid)
    # No SNP information, add as NA
    metho$issnp = data.table(x=rep(NA, nrow(metho$pos)))
    names(metho$issnp) = sid
    return(metho)
}

read_normal_as_metho <- function(sample_id, ac_file){
    ac = data.table::fread(ac_file, verbose=F, showProgress=F)
    ac = ac[!is.na(m),]

    pos = ac[, .(chrom, start, end)]
    meth = ac[, .(m)]
    cov = ac[, .(total_counts_m)]
    issnp = ac[, .(is_snp=!is.na(POS))]

    # Name after sample
    names(meth) = sample_id
    names(cov) = sample_id
    names(issnp) = sample_id

    metho = list(pos=pos, meth=meth, cov=cov, issnp=issnp)

    # Drop any CpGs with duplicated positions in manifest (old allele counter)
    keep_pos = !duplicated(pos)
    metho = lapply(metho, function(x) x[keep_pos,])
    return(metho)
}

# Merge SINET methylation objects
merge_metho <- function(x, y) {
    # x and y must have pos/meth/cov data tables in a list
    # ensure pos are characters
    x$pos$chrom = as.character(x$pos$chrom)
    y$pos$chrom = as.character(y$pos$chrom)

    # Format for overlap
    xp <- x$pos
    xp$xindex <- seq_len(nrow(xp))
    yp <- y$pos
    yp$yindex <- seq_len(nrow(yp))
    setkey(xp, chrom, start, end)
    setkey(yp, chrom, start, end)

    # Set new CpG positions
    pos <- merge(xp, yp, all = T)
    # As some keys are duplicated in some CAMDAC files, we remove them here
    pos <- pos[!duplicated(pos, by = c("chrom", "start", "end"))]

    # Set feature vectors for new methylation data
    featvecs <- setdiff(names(x), "pos")
    obj = list("pos"=pos)
    for(feat in featvecs){
        xfeat = x[[feat]]
        yfeat = y[[feat]]
        obj[[feat]] = cbind(
            xfeat[pos$xindex, ],
            yfeat[pos$yindex, ]
        )
        names(obj[[feat]]) = c(names(xfeat), names(yfeat))
    }

    # Format obj pos
    obj$pos = obj$pos[, .(chrom, start, end)]

    return(obj)
}


metho_make <- function(chrom, start, end, meth, cov, sname){
    obj <- list(
        pos = data.table(chrom, start, end),
        meth = data.table(x=meth),
        cov = data.table(x=cov)
    )
    names(obj$meth) <- sname
    names(obj$cov) <- sname
    return(obj)
}

read_pure_metho <- function(sample_id, pure_file){
    ac = data.table::fread(pure_file, verbose=F, showProgress=F)
    ac = ac[!is.na(m) | !is.na(m_t),]

    pos = ac[, .(chrom, start, end)]

    meth = ac[, .(m)]
    meth_pure = ac[, .(m_t)]

    cov = ac[, .(cov)]
    cov_pure = ac[, .(cov_t)]

    issnp = ac[, .(is_snp=SNP==1)]

    # Name after sample
    names(meth) = sample_id
    names(meth_pure) = sample_id
    names(cov) = sample_id
    names(cov_pure) = sample_id
    names(issnp) = sample_id

    metho = list(pos=pos, meth=meth, meth_pure=meth_pure, cov=cov, cov_pure=cov_pure, issnp=issnp)

    # Drop any CpGs with duplicated positions in manifest (old allele counter)
    keep_pos = !duplicated(pos)
    metho = lapply(metho, function(x) x[keep_pos,])
    return(metho)
}