# This file contains scripts for validating panels as proxy for adjacent normals
renv::load("~/projects/wgbs_camdac_dev") # Skip local renv to use CAMDAC DMP function
devtools::load_all("~/projects/VanLoo_CAMDAC")
library(data.table)
library(dplyr)
library(logger)
library(foreach)
library(doParallel)
library(future)
library(listenv) # Req for future
library(matrixStats)
library(future.batchtools)

dispatch_panel_pipeline = function(infile, exp){
    # Load methylation object
    meth_obj = qs::qread(infile)


    if(grepl('validation', exp)){
        val_obj = pipe_panel_validation(meth_obj)
        return(val_obj)
    }else if(grepl('saturation', exp)){
        # Downsampled to 2.5m CpGs for speed
        meth_obj = lapply(meth_obj, function(x) x[1:2500000,])
        sat_obj = pipe_panel_saturation(meth_obj)
        return(sat_obj)
    }

    return(NA)
}

# Create a pipeline that calculates stats for each sample relative to a panel of others
#' @param pandata A dataframe panel list of normal methylomes. See `dummy_panel_data()` for an example.
#  @param samples A vector of samples to use for analysis. Can be NULL.
pipe_panel_validation <- function(pandata, samples = NULL) {

    # Subset samples if requested
    if (!is.null(samples)) {
        pandata <- subset_samples(pandata, samples)
    }

    # Validate panel data
    validate_pandata(pandata)

    # Get list of samples to work over
    csamples <- list_panel_samples(pandata)

    # Iterate over each sample, creating a panel and running the pipeline
    presult <- list()
    for (item in csamples) {
        log_info("Running pipeline for {item}")

        # Create panel from all other samples
        pother <- setdiff(csamples, item)
        pother_name <- paste0(item, "_panel")
        panel <- make_panel(pandata$meth, pandata$cov, pother, pother_name)

        # Get panel DMP stats with normal
        pdmo <- call_dm_obj(pandata, item, panel, pother_name)
        pdm <- call_dm_stats(pdmo)

        # Add annotations
        pdm$item <- item
        #pdm$pother <- paste0(pother, collapse = "|") # Too much data. Rerun if provenance needed as output file.
        itemix <- which(item == list_panel_samples(pandata))
        pdm$tissue <- pandata$tissue[itemix]
        pdm$sex <- pandata$sex[itemix]
        pdm$assay <- pandata$assay[itemix]
        pdm$source <- pandata$source[itemix]

        # Format pdmo
        non_chr_names = setdiff(names(pdmo), c("chrom", "start", "end"))
        pdmo = data.table(pdmo)
        pdmo = pdmo[!duplicated(pdmo[, .(chrom,start,end)])]
        for(i in non_chr_names){
            setnames(pdmo, i, paste0(item, "..", i))
        }

        presult[[item]] <- list(pdm=pdm, pos=pdmo)
    }

    # Merge the separate results into objects and return
    result <- data.table::rbindlist(lapply(presult, function(x) x$pdm))
    result_pos <- Reduce(function(x, y) merge(x, y, by=c("chrom", "start", "end"), all=T), lapply(presult, function(x) x$pos))

    comp = list(result=result, pos=result_pos)

    return(comp)
}


pipe_panel_saturation <- function(pandata, nsubset=50){
    # Saturation analysis to determine panel performance using N samples to create the panel

    N = c( seq(10), seq(15,30,5), c(40, 50, 75, 100)) # Number of items in the panel
    csamples <- list_panel_samples(pandata)
    if(!is.null(nsubset)){
        csamples = sample(csamples, min(nsubset, length(csamples)-1))
    }

    # Function that takes a single N and runs the panel pipeline for K panels
    all_res = list()
    for (current_n in N) {
        # N.B. we don't include panel of all other samples as calculated in the validation pipeline
        if(current_n >= length(csamples)-1){
            log_info("Skipping N={current_n} as there are not enough samples")
            next
        }

        # Run for this N
        doParallel::registerDoParallel(cores = 15)
        N_res = foreach(item=csamples, .combine='rbind') %dopar% {
            log_info("Running SATURATION pipeline for {item}")

            # Create as many panel(s) of size N as possible
            pother <- setdiff(csamples, item)
            pother_panels <- make_n_panels(pandata$meth, pandata$cov, pother, current_n, maxK=5)
            pother_name <- paste0(item, "_panel")

            # Get panel DMP stats with normal
            pdm_res <- lapply(
                seq_along(pother_panels),
                function(pp){
                    panel_ix = names(pother_panels)[pp]
                    panel_obj = pother_panels[[panel_ix]]
                    panel_name = list_panel_samples(panel_obj)
                    pdm_n = call_dm_obj(pandata, item, panel_obj, panel_name)
                    pdm_n = call_dm_stats(pdm_n)
                    pdm_n$panel_ix = panel_ix
                    return(pdm_n)
                }
            )
            pdm_res = data.table::rbindlist(pdm_res)
            pdm_res$item = item
            pdm_res$current_n = current_n

            # Add to cache
            return(pdm_res)
        }
        doParallel::stopImplicitCluster()

        all_res[[as.character(current_n)]] = N_res
    }
    psat = data.table::rbindlist(all_res)
    return(psat)


}

# Create a dummy panel data object for testing with 5 samples and 1000 CpGs
dummy_panel_data <- function(samps=5, n_cpg=1000) {
    pandata <- list(
        pos = data.frame(
            chrom = rep("chr1", n_cpg),
            start = seq(1, n_cpg),
            end = seq(1, n_cpg) + 1
        ),
        meth = matrix(rbeta(n_cpg * samps, 1, 1), nrow = n_cpg, ncol = samps) %>%
            as.data.frame() %>%
            setnames(paste0("P", 1:samps)),
        cov = matrix(rbinom(n_cpg * samps, 100, 0.5), nrow = n_cpg, ncol = samps) %>%
            as.data.frame() %>%
            setnames(paste0("P", 1:samps))#,
        # matched = sample(c(T, T, T, F, T), samps, replace = T),
        # tissue = rep("prostate", samps),
        # sex = rep("XY", samps),
        # assay = rep("WGBS", samps),
        # source = rep("SU2C", samps)
    )
    return(pandata)
}

# Validate panel data object. Namely check that:
# - names of meth and cov are equivalent
#  - number of columns in meth is the same as length of matched, tissue, sex, assay and source
validate_pandata <- function(pandata) {
    if (!all(names(pandata$meth) == names(pandata$cov))) {
        stop("Names of meth and cov are not equivalent")
    }
    # if (ncol(pandata$meth) != length(pandata$matched) ||
    #     ncol(pandata$meth) != length(pandata$tissue) ||
    #     ncol(pandata$meth) != length(pandata$sex) ||
    #     ncol(pandata$meth) != length(pandata$assay) ||
    #     ncol(pandata$meth) != length(pandata$source)) {
    #     stop("Number of columns in meth is not the same as length of matched, tissue, sex, assay and source")
    # }
    return(TRUE)
}

list_panel_samples <- function(x) names(x$meth)

make_panel <- function(meth, cov, samples, newname) {
    # Create a panel by merging data from the and cov dataframes
    # Assumes the same CpGs are in the same location for each row

    meths <- dplyr::select(meth, all_of(samples))
    covs <- dplyr::select(cov, all_of(samples))


    panel <- list(
        meth = data.table(x = rowMedians(as.matrix(meths), na.rm = T)),
        cov = round(data.table(x = rowMedians(as.matrix(covs), na.rm = T)),0)
    )
    names(panel$meth) <- newname
    names(panel$cov) <- newname

    return(panel)
}

call_dm_obj <- function(obj1, sample1, obj2, sample2){
    # Setup DMP data
    # mbdiff is theoretically tumor - normal, making tumor=sample1 and normal=sample2
    mbdiff <- obj1$meth[[sample1]] - obj2$meth[[sample2]]

    M_n <- round(obj2$meth[[sample2]] * obj2$cov[[sample2]], 0)
    UM_n <- obj2$cov[[sample2]] - M_n
    M <- round(obj1$meth[[sample1]] * obj1$cov[[sample1]], 0)
    UM <- obj1$cov[[sample1]] - M

    # Call dmps
    log_info("Calling DMPs for {sample1} vs {sample2}")
    pdmp <- dmp_call_pipe(mbdiff, M_n, UM_n, M, UM,
        mtdiff = NULL,
        effect_size = 0.2, prob = 0.99, itersplit = 1e6, ncores=1
    )
    pdmp <- cbind(obj1$pos, pdmp)
    pdmp <- pdmp[!is.na(pdmp$mbdiff), ]
    pdmp = pdmp %>%
        dplyr::mutate(DMP_b = ifelse(is.na(DMP_b),
            ifelse(abs(mbdiff) > 0.2, "N.S.", NA_character_), DMP_b))

    return(pdmp)
}

call_dm_stats <- function(pdmp) {

    # Flatten to stats
    statdmp <- pdmp  %>%
        dplyr::count(DMP_b) %>%
        dplyr::mutate(categories = dplyr::case_when(
            is.na(DMP_b) ~ "non",
            DMP_b=="N.S." ~ "ns",
            TRUE ~ DMP_b )) %>%
        dplyr::mutate(categories = paste0("perc_", categories)) %>%
        dplyr::select(-DMP_b) %>%
        dplyr::mutate(n_cg = sum(n)) %>%
        dplyr::mutate(percs = round((n * 100) / n_cg, 2)) %>%
        tidyr::pivot_wider(
            id_cols = c("n_cg"),
            names_from = "categories", values_from = "percs"
        )

    # Calculate additional stats
    statdmp$mdiff_mean <- mean(pdmp$mbdiff, na.rm = T)
    statdmp$mdiff_iqr <- IQR(pdmp$mbdiff, na.rm = T)
    statdmp$euclidean <- sqrt(sum(pdmp$mbdiff^2, na.rm = T))

    return(statdmp)
}

subset_samples <- function(pandata, samples) {
    curr_samples <- list_panel_samples(pandata)
    sel_ix <- match(samples, curr_samples)
    pandata <- lapply(
        pandata,
        function(x) {
            if (is.data.frame(x)) {
                x <- x[, sel_ix]
            } else {
                x <- x[sel_ix]
            }
            return(x)
        }
    )
    return(pandata)
}


# Returns a list of panels made up of current_n samples
make_n_panels <- function(meth, cov, pother, current_n, maxK=5){
    # Calculate at least 5 if possible

    # Get combinations of samples at current_n
    # If too many combinations (> 5e6) then simply sample
    ncombinations = choose(length(pother), current_n)
    if(ncombinations > 1e5){
        log_warn("Combinations for panel size {current_n} in data {length(pother)} are > 1e5. Using sampling protocol.")
        col_sel = seq(min(current_n, maxK))
        combos = lapply(col_sel, function(i) sample(pother, current_n))
    } else {
        combos = combn(pother, current_n)
        ncombos = ncol(combos)
        col_sel = sample(seq(ncombos), min(maxK, ncombos))
        combos = lapply(col_sel, function(i) combos[,i])
    }

    # Create a panel for each combination
    mpanels = lapply(
        seq_along(combos),
        function(i){
            nother = combos[[i]]
            nother_name = i
            panel = make_panel(meth, cov, nother, as.character(nother_name))
        }
    )
    names(mpanels) = col_sel
    return(mpanels)
}
