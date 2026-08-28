# Figure S5F
## Aim: Use HDI estimates from simulations (drawn from RRBS parameters) to show how CAMDAC 
##    deconvolution power varies by purity, copy number and coverage

# Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
source('src/power_estimates.R')
options(bitmapType='cairo')

# Load initial tracerX simulation data
# Source: /nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/ECADIEUX_FULL/txMeth/power_sim_results_list.RData
load('data/pre/simulations/power_sim_results_list.RData')
sim = results_list

# Convert HDI estimates (? per simulated CpG allele ?) into a long-format clean data.table
hdi = data.table::rbindlist(
    lapply(
        seq_along(sim),
        function(i) parse_sim_hdi_long(sim[[i]], i)
    )
)

# Cleanup HDI table with a number of filters:
# - Remove rows with missing mt_hdi_width values
hdi = hdi[!is.na(mt_hdi_width)]
# - Optional: Keep only non-DMP (best-estimate HDI)
hdi = hdi[DMP_true=='no']

# Setup datasets to observe relationship between HDI width and parameters
#   Note: Datasets designed to 'fix' the parameter of interest
d_coverage = hdi[CN==2 & dplyr::between(purity, 0.59-0.1, 0.59+0.1)] # +/- 10% of median purity (0.59)
d_purity = hdi[CN==2 & dplyr::between(mt_cov, 67-10, 67+10),] # +/- 10X of median coverage (67X)
d_CN = hdi[dplyr::between(purity, 0.59-0.2, 0.59+0.2) & dplyr::between(mt_cov, 67-20, 67+20),] # +/- 10% of median purity (0.59) and +/- 10X of median coverage (67X)

# Additional cleanup
## Cap coverage and select non-DMPs
d_coverage = d_coverage[mt_cov <= 300]
d_coverage = d_coverage[DMP_true=='no'] # Observed DMP sites have higher HDI
## Set purity as a categorical value with cut
# d_purity$p_cut = cut(
#     d_purity$purity,
#     breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1),
#     labels=c('0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-1')
# )
d_purity$p_cut = cut(
    d_purity$purity,
    breaks=c(0, 0.25, 0.5, 0.75, 1),
    labels=c('0-0.25', '>0.25-0.5', '>0.5-0.75', '>0.75-1')
)
## When copy number values are above 6, set as '6+' and limit to below 10
d_CN = d_CN[d_CN$CN < 10,]
d_CN$CN = ifelse(d_CN$CN >= 6, '6+', as.character(d_CN$CN))


# Generate plots against the HDI parameter for each of these datasets
expr = paste0("HDI^99 range")
p_coverage = (
    ggplot(d_coverage, aes(x=mt_cov, y=mt_hdi_width)) +
    geom_point(alpha=0.05, fill='lightgrey', shape=21, size=0.2) +
    geom_smooth(method='loess', color='black') +
    labs(title='', x='Coverage', y=expression(paste(HDI^99, " range")))    +
    theme_minimal() +
    scale_y_continuous(limits=c(0, 1))
)

p_purity = (
    ggplot(d_purity, aes(x=p_cut, y=mt_hdi_width)) +
    geom_jitter(alpha=0.05, size=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.9) +
    labs(title='', x='Purity', y=expression(paste(HDI^99, " range"))) +
    theme_minimal() +
    scale_y_continuous(limits=c(0, 1)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) 
    # theme(axis.text.y=element_blank()) +
    # labs(y="")
)

p_CN = (
    ggplot(d_CN, aes(x=CN, y=mt_hdi_width)) +
    geom_jitter(alpha=0.05, size=0.1) +
    geom_boxplot(outlier.shape=NA, alpha=0.9) +
    geom_smooth(method='lm') +
    labs(title='', x='Copy number', y=expression(paste(HDI^99, " range"))) +
    scale_y_continuous(limits=c(0, 1)) +
    theme_minimal()
    # theme(axis.text.y=element_blank()) +
    # labs(y="")
)

# Save plots
p_patch = patchwork::wrap_plots(p_coverage, p_purity, p_CN, nrow=1)
#ggsave('figures/power_analysis_hdi_simulations.png', p_patch, width=10, height=5, units='in', dpi=300)

ggsave('figures/fig_s5_power_analysis_hdi_simulations.png', p_patch, width=8.27, height=3.69, units='in', dpi=400)



