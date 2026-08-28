# Libraries ----
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(qs)
library(fs)
library(logger)
options(bitmapType = "cairo")

#' Pipeline to measure confounding in bulk tumors
#' pum: CAMDAC pure meth object at sites unmethylated in normal
#' phm: CAMDAC pure meth object at sites hypermethylated in normal
#' pid: patient ID
#' csvdir: directory to save CSVs
#' figdir: directory to save figures
pipeline_confounding <- function(pum, phm, pid, csvdir, figdir) {
    log_info("Running confounding pipeline for {pid}")
    csvhm <- fs::path(csvdir, paste0(pid, "_HM.csv"))
    csvnu <- fs::path(csvdir, paste0(pid, "_NU.csv"))

    if (all(fs::file_exists(c(csvhm, csvnu)))) {
        log_info("CSVs already exist. Skipping.")
        return()
    }

    log_info("Creating histograms")
    # Make histogram data tables
    dt_nu <- make_hist_dt_cn(pum, id = paste0(pid, "_NU")) # normal unmeth
    dt_hm <- make_hist_dt_cn(phm, id = paste0(pid, "_HM")) # normal hypermeth

    log_info("Plotting histograms by CN")
    # Plot CN states and peak solution
    p_nu <- plot_histdt(dt_nu, peak = "max")
    p_hm <- plot_histdt(dt_hm, peak = "min")
    p_combined <- p_nu / p_hm
    figfile <- fs::path(figdir, paste0(pid, ".png"))
    ggsave(figfile, p_combined, width = p_nu$pwidth*1.5, height = p_nu$pheight + p_hm$pheight)

    log_info("Saving solutions to CSV")
    # Save solutions to a CSV file
    csv_nu <- histdt_to_csv(dt_nu, type = "max", outfile = csvnu)
    csv_hm <- histdt_to_csv(dt_hm, type = "min", outfile = csvhm)

    log_info("Complete:: {pid}.")
}

gather_confounding_csvs <- function(csvdir) {
    # Plot confounding summary
    conf_csvs <- fs::dir_ls(csvdir, glob = "*.csv")
    dconf <- rbindlist(lapply(conf_csvs, fread))
    dconf <- dconf[qc==0,]
    dconf
}


# Function to calculate histogram for tumor sample
#' m: methylation vector
#' breaks: breaks for histogram
calc_meth_hist <- function(m, breaks = NULL) {
    if (is.null(breaks)) {
        breaks <- seq(from = 0, to = 1, by = 0.05)
    }
    histdata <- hist(m, breaks = breaks, plot = FALSE)

    # Cap histogram at upper quantile to avoid outlies
    uq <- quantile(histdata$counts, 0.9)
    histdata$cap <- ifelse(histdata$counts > uq, uq, histdata$counts)

    return(histdata)
}

# Function to save data as data.table
#' e: output from calc_meth_hist
meth_hist_to_dt <- function(e) {
    data.table(
        field = e$xname,
        breaks = list(e$breaks),
        counts = list(e$counts),
        cap = list(e$cap),
        density = list(e$density)
    )
}


# Function to find peaks
# source: https://github.com/stas-g/findPeaks/tree/master
#' x: numeric vector of counts
#' m: number of points on each side of the peak
find_peaks <- function(x, m = 3) {
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i) {
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) {
            return(i + 1)
        } else {
            return(numeric(0))
        }
    })
    pks <- unlist(pks)
    pks
}


# Pipeline from methylation to histogram data
#' x: data.table with columns CN, p, m
#' id: id for sample
make_hist_dt_cn <- function(x, id = "") {
    # Check input
    stopifnot(all(c("CN", "p", "m") %in% names(x)))

    # Select CN states with at least 1000 representative CpGs for analysiss
    CN_states <- unique(x$CN)
    CN_states <- CN_states[sapply(CN_states, function(i) sum(x$CN == i) > 1000)]

    cnhists <- lapply(CN_states, function(i) {
        # Get histogram
        hhist <- calc_meth_hist(x[CN == i, ]$m)
        # Convert to data.table
        dhist <- meth_hist_to_dt(hhist)
        #  Find peaks
        dhist$peakix <- list(find_peaks(hhist$cap, m = 3))
        dhist$peakbin <- list(hhist$breaks[dhist$peakix[[1]]])

        # Add max and min peaks
        dhist$maxpeak <- max(dhist$peakbin[[1]])
        dhist$minpeak <- min(dhist$peakbin[[1]])

        # Add CN state, purity and id for posterity
        dhist$CN <- i
        dhist$CN_state <- ifelse(dhist$CN>=5, "5+", as.character(dhist$CN))
        dhist$id <- id
        dhist$purity <- x$p[[1]]
        return(dhist)
    })
    return(data.table::rbindlist(cnhists))
}


# Function to plot the histograms from each row of histdt
#' histdt: output from make_hist_dt_cn
#' peak: which peak to plot (max or min)
plot_histdt <- function(histdt, peak = "max") {
    # For each row, plots capped methylation as histogram
    # Adds blue dashed line for purity
    # Adds red solid line for current 'peak' solution
    nplots <- nrow(histdt)
    peak_plots <- lapply(seq(nplots), function(i) {
        plot_histdt_row(histdt[i, ], peak = peak)
    })

    ncol <- 6
    nrow <- ceiling(nplots / ncol)
    if (nrow == 1) {
        ncol <- nplots
    }
    peak_patch <- patchwork::wrap_plots(peak_plots, ncol = ncol)
    peak_patch$pwidth <- 3 * ncol
    peak_patch$pheight <- 3 * nrow

    return(peak_patch)
}

#' Function to plot a single row of histdt
#' row: a single row of histdt
#' peak: which peak to plot (max or min)
plot_histdt_row <- function(row, peak = "max") {
    stopifnot(peak %in% c("max", "min"))
    counts <- unlist(row$cap)
    breaks <- unlist(row$breaks)
    breaks <- breaks[1:length(breaks) - 1]
    purity <- ifelse(peak == "max", row$purity, 1 - row$purity)
    peakid <- paste0(peak, "peak")
    ptitle <- paste0(row$id, "::CN=", row$CN)

    ggplot(data.frame(x = breaks, y = counts), aes(x = x, y = y)) +
        geom_col() +
        geom_vline(xintercept = purity, color = "blue", linetype = "dashed") +
        geom_vline(xintercept = row[[peakid]], col = "red") +
        labs(title = ptitle)
}


# Function to write a CSV with solutions
#' x: output from make_hist_dt_cn
#' outfile: path to output csv file
histdt_to_csv <- function(x, type, outfile = NULL) {
    id <- x$id
    purity <- x$purity
    CN <- x$CN
    CN_state <- x$CN_state
    nfield <- x$field # The methylation field used ('m')
    nbreaks <- sapply(x$breaks, function(i) i[2] - i[1]) # Get the bin width used
    peakbins <- sapply(x$peakbin, function(i) paste0(i, collapse = ";"))
    peaktype <- type
    maxpeak <- x$maxpeak
    minpeak <- x$minpeak
    qc <- 0

    # Create data table from the above
    dt <- data.table(
        id = id,
        purity = purity,
        CN = CN,
        CN_state = CN_state,
        nfield = nfield,
        nbreaks = nbreaks,
        peakbins = peakbins,
        peaktype = peaktype,
        maxpeak = maxpeak,
        minpeak = minpeak,
        qc = qc
    )

    # Write output file
    if (!is.null(outfile)) {
        data.table::fwrite(dt, outfile)
    }

    return(dt)
}

# Function to write create the confounding summary plot given a CSV with multiple solutions
#' x: output from histdt_to_csv
plot_conf_cohort_peak <- function(xin, peaktype = "max") {
    # Set the type of plot we are creating (sites at unmeth normal [max] or meth normal [min])
    if (peaktype == "max") {
        x <- xin[xin$peaktype == "max", ]
        # xlab <- "WGS-derived tumor purity"
        # xval <- x$purity
        yval <- "maxpeak"
        ylab <- "Hypermethylation peak at CpGs\n unmethylated in normal"
    } else if (peaktype == "min") {
        x <- xin[xin$peaktype == "min", ]
        # xlab <- "1 - (WGS-derived tumor purity)"
        # xval <- 1 - x$purity
        yval <- "minpeak"
        ylab <- "Hypomethylation peak at CpGs\n methylated in normal"
    } else {
        stop("peaktype must be max or min")
    }

    # Set plot
    ggplot(x, aes(x = purity, y = get(yval))) +
        geom_point(aes(color=factor(CN_state))) +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 1)) +
        geom_smooth(aes(color = factor(CN_state)), method = lm, se = FALSE) +
        scale_color_brewer(palette = "Set2") +
        theme_classic() +
        labs(color = "Total CN", y = ylab, x = "Tumor purity")
}

# Function to extract identifier from pure methylation file
extract_pure_meth_fid <- function(x) {
    basename <- fs::path_file(x)
    str_remove(basename, "\\.pure.*")
}


plot_conf_dot <- function(clnconf, conf_dot_plot){

    # Plot confounding summary over cohort
    maxconf <- plot_conf_cohort_peak(clnconf, peaktype = "max")
    minconf <- plot_conf_cohort_peak(clnconf, peaktype = "min")
    pconf = maxconf / minconf
    ggsave(conf_dot_plot, pconf, width = 5, height = 8)
    return(pconf)

}

plot_conf_tile <- function(clnconf, conf_tile_plot){
    # Plot confounding summary over cohort
    maxconf <- plot_conf_cohort_tile(clnconf, peaktype = "max")
    minconf <- plot_conf_cohort_tile(clnconf, peaktype = "min")
    pconf = maxconf / minconf
    ggsave(conf_tile_plot, pconf, width = 7, height = 8)
    return(pconf)

}

plot_conf_cohort_tile <- function(clnconf, peaktype){
    stopifnot(peaktype %in% c("max", "min"))
    peakname = paste0(peaktype, "peak")

    gradcolor = ifelse(peaktype == "max", "darkred", "white")
    gradlow = ifelse(peaktype == "max", "white", "aquamarine4")
    # Select give peak type
    sel = clnconf$peaktype == peaktype
    x = clnconf[sel,]

    # Create tile plot
    x$id = forcats::fct_reorder(x$id, x$purity)
    ptile = ggplot(x, aes(x=id, y=CN_state, fill=get(peakname))) + 
        geom_tile() + scale_fill_gradient(low=gradlow, high=gradcolor) + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, size=4)) +
        labs(fill=peakname)
    pbar = ggplot(unique(x[, .(id, purity)]), aes(x=id, y=purity)) + geom_col() + theme_classic() +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    pcomb = patchwork::wrap_plots(pbar, ptile, heights=c(0.1,1))
    return(pcomb)
}