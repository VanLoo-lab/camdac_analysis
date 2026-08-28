# Aim: Plot p values using purity / parallel
library(data.table)
library(ggplot2)
library(tidyr)
library(stringr)
library(fs)
library(qs)
library(logger)
library(foreach)
library(doParallel)
library(patchwork)
library(figpatch)
options(bitmapType = "cairo")

# Functions ----

parse_ora_f <- function(f) {
    e <- fs::path_file(f) %>% stringr::str_remove(".qs")
    e <- stringr::str_split(e, "_") %>% unlist()
    dt <- qs::qread(f)$ora
    dt$source <- e[[2]]
    dt$purityf <- e[[3]]
    return(dt)
}

get_ora_data <- function(outdir, suffix, src) {
    ora_files <- fs::dir_ls(outdir, glob = paste0("*ORA_*", suffix, ".qs"))
    ores <- data.table::rbindlist(lapply(ora_files, parse_ora_f))
    ores <- dplyr::filter(ores, source == src)
    ores <- clean_ores(ores)
    return(ores)
}

clean_ores <- function(ores) {
    ores %>%
        dplyr::mutate(lpadj = -log10(padj)) %>%
        dplyr::mutate(index = 1:nrow(.)) %>%
        dplyr::select(pathway, dtype, hallmark, source, dset, lpadj) %>%
        tidyr::pivot_wider(names_from = dset, values_from = lpadj) %>%
        dplyr::mutate(dtype = ifelse(dtype == "hyper", "Hypermethylation", "Hypomethylation"))
}

plot_ora_data <- function(ores) {
    # Set equal limits
    mylim = max(c(ores$bulk, ores$pure))
    ggplot(ores, aes(x = bulk, y = pure, fill = hallmark)) +
        geom_point(shape = 21, alpha = .5) +
        facet_wrap(~dtype, ncol = 1, scales = "free") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        theme_bw() +
        theme(strip.background = element_blank()) +
        labs(x = "Bulk tumour\n-log10[q-value]", y = "CAMDAC-purified\n-log10[q-values]", fill = "MSigDb Collection") +
        scale_y_continuous(n.breaks=6) +
        scale_fill_manual(labels = c("hall" = "Hallmark Pathways", "o3ca" = "Curated Cancer\nCell Atlas", "onc" = "Oncogenic\nSignature Genes"), values = c("purple", "black", "orange"))
}

# Main ----
outdir <- fs::dir_create("data/pre/gsea_500") # Source data dir, not true outdir
figdir <- fs::dir_create("figures/241219_ora")

# Create runsheet
suffix <- c("lq", "hq") # lower quantile and upper quantile based on median
source <- c("luad", "lusc", "wgbs") # datasets
rs <- expand.grid(suffix = suffix, source = source)

# Generate plot for each record
registerDoParallel(cores = 4)
ora_plot_parsed <- foreach(i=1:nrow(rs)) %dopar% {

    # Setup record
    suffix <- rs$suffix[i]
    src <- rs$source[i]
    ora_id = paste0(src, "_", suffix)
    log_info("Generating ORA plots for {ora_id}")

    # Get data and cache
    ores <- get_ora_data(outdir, suffix, src)

    # Get plot and cache for hyper
    dty="Hypermethylation"
    p_oreshy <- plot_ora_data(ores %>% dplyr::filter(dtype==dty & hallmark=='onc'))
    # Save singular plot 
    ggsave(fs::path(figdir, paste0("osingle_", src, "_", suffix, '_', dty, ".png")), p_oreshy, width = 3, height = 3, dpi = 300)
    ggsave(fs::path(figdir, paste0("osingle_", src, "_", suffix, '_', dty, ".pdf")), p_oreshy, width = 3, height = 3, dpi = 300)

    # Get plot and cache for hypo
    dty="Hypomethylation"
    p_oresho <- plot_ora_data(ores %>% dplyr::filter(dtype==dty & hallmark=='onc'))
    # Save singular plot 
    ggsave(fs::path(figdir, paste0("osingle_", src, "_", suffix, '_', dty, ".png")), p_oresho, width = 3, height = 3, dpi = 300)
    ggsave(fs::path(figdir, paste0("osingle_", src, "_", suffix, '_', dty, ".pdf")), p_oresho, width = 3, height = 3, dpi = 300)

    # Cache output
    res <- list(
        ora_id=ora_id,
        ores=ores,
        p_oreshy=p_oreshy,
        p_oresho=p_oresho
    )

    return(res)
}
names(ora_plot_parsed) <- sapply(ora_plot_parsed, function(x) x$ora_id)
qs::qsave(ora_plot_parsed, fs::path(figdir, "ora_plot_parsed.qs"))

# Load plots
orp <- qs::qread(fs::path(figdir, "ora_plot_parsed.qs"))
hy_lq_plots = list(
    orp[["luad_lq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUAD"),
    orp[["lusc_lq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUSC", y=""),
    orp[["wgbs_lq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "mCRPC", y="")
)
ho_lq_plots = list(
    orp[["luad_lq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUAD"),
    orp[["lusc_lq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUSC", y=""),
    orp[["wgbs_lq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "mCRPC", y="")
)
hy_hq_plots = list(
    orp[["luad_hq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUAD"),
    orp[["lusc_hq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUSC", y=""),
    orp[["wgbs_hq"]]$p_oreshy + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "mCRPC", y="")
)
ho_hq_plots = list(
    orp[["luad_hq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUAD"),
    orp[["lusc_hq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "LUSC", y=""),
    orp[["wgbs_hq"]]$p_oresho + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = "mCRPC", y="")
)

# Setup patches
hy_lq_patch = patchwork::wrap_plots(hy_lq_plots, ncol = 3)
ho_lq_patch = patchwork::wrap_plots(ho_lq_plots, ncol = 3)
hy_hq_patch = patchwork::wrap_plots(hy_hq_plots, ncol = 3)
ho_hq_patch = patchwork::wrap_plots(ho_hq_plots, ncol = 3)

# Plot pngs
ggsave(fs::path(figdir, "ora_hy_lq_patch.png"), hy_lq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_ho_lq_patch.png"), ho_lq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_hy_hq_patch.png"), hy_hq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_ho_hq_patch.png"), ho_hq_patch, width = 5.5, height = 3, dpi = 300)

# Plot pdfs
ggsave(fs::path(figdir, "ora_hy_lq_patch.pdf"), hy_lq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_ho_lq_patch.pdf"), ho_lq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_hy_hq_patch.pdf"), hy_hq_patch, width = 5.5, height = 3, dpi = 300)
ggsave(fs::path(figdir, "ora_ho_hq_patch.pdf"), ho_hq_patch, width = 5.5, height = 3, dpi = 300)
