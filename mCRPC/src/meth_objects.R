library(data.table)
# Helper function for SINET methylation objects

metho_make <- function(chrom, start, end, meth, cov, sname){
    obj <- list(
        cg = data.table(chrom, start, end),
        meth = data.table(x=meth),
        cov = data.table(x=cov)
    )
    names(obj$meth) <- sname
    names(obj$cov) <- sname
    return(obj)
}

# Merge SINET methylation objects
metho_merge <- function(x, y) {
    # x and y must have pos/meth/cov data tables in a list

    # Format for overlap
    xp <- x$cg
    xp$xindex <- seq_len(nrow(xp))
    yp <- y$cg
    yp$yindex <- seq_len(nrow(yp))
    setkey(xp, chrom, start, end)
    setkey(yp, chrom, start, end)

    # Set new CpG positions
    pos <- merge(xp, yp, all = T)
    # As some keys are duplicated in some CAMDAC files, we remove them here
    pos <- pos[!duplicated(pos, by = c("chrom", "start", "end"))]

    # Set new methylation data
    meth <- cbind(
        x$meth[pos$xindex, ],
        y$meth[pos$yindex, ]
    )

    # Set new coverage data
    cov <- cbind(
        x$cov[pos$xindex, ],
        y$cov[pos$yindex, ]
    )

    # Set cg sites
    cg <- pos[, .(chrom, start, end)]

    # Build output object
    obj <- list(cg = cg, meth = meth, cov = cov)
    return(obj)
}

metho_filter_samples <-function(meth, sid) {
    insamples = names(meth$meth)
    stopifnot(all(sid %in% insamples))
    meth$cov = data.table(meth$cov)[, ..sid]
    meth$meth = data.table(meth$meth)[, ..sid]
    meth$cg = data.table(meth$cg)
    return(meth)
}

metho_rename_sample <- function(meth, sid, newname){
    insamples = names(meth$meth)
    stopifnot(sid %in% insamples)
    nameix = which(insamples==sid)
    names(meth$meth)[nameix] = newname
    names(meth$cov)[nameix] = newname
    return(meth)
}

metho_filter_chrom <- function(x, cfilt){
    sel = x$cg$chrom == cfilt
    mf = lapply(x, function(y) y[sel, ])
    return(mf)
}

metho_filter_region <- function(x, chrom, start, end){
    sel = x$cg$chrom == chrom & x$cg$start >= start & x$cg$end <= end
    mf = lapply(x, function(y) y[sel, ])
    return(mf)
}

metho_filter_bool <- function(x, bool){
    sel = bool
    mf = lapply(x, function(y) y[sel, ])
    return(mf)
}

metho_anno_cpgs <- function(x, annos){
    # Annos should be a data.table with chrom, start, end.
    # Returns object with x$cg annotated
    cg = x$cg
    setkey(cg, chrom, start, end)
    setkey(annos, chrom, start, end)
    
    cgo = foverlaps(cg, annos, nomatch=NA, mult="first")
    cgo[, cgi_start:=start]
    cgo[, cgi_end:=end]
    cgo[, start:=i.start]
    cgo[, end:=i.end]
    cgo[, `:=`(i.start=NULL, i.end=NULL)]
    
    stopifnot(nrow(x$cg) == nrow(cgo))
    x$cg = cgo
    return(x)
}

metho_select_samples <- function(x, snames){
    x$meth = x$meth[, ..snames]
    x$cov = x$cov[, ..snames]
    return(x)
}

metho_sample_n <- function(x, n){
    recs = nrow(x$meth)
    sampe_vec <- sample(1:recs, min(n, recs))
    x2 <- lapply(
        x, 
        function(y) y[sampe_vec, ]
    )
    return(x2)
}

metho_count_cg <- function(x){
    ncg = data.table(
        sample = names(x$cov),
        n_cpg = colSums(!is.na(x$cov))
    )
}