# Reviewer comments require us to reanalyse the FACS tumor data
## 1C9 questions the 77% agreement between CAMDAC-mt DMPs and FACS-mt DMPs.
## Is this a good agreement? We therefore demonstrate that this is an improvement over the bulk.

## Lili performed this comparison it seems.
fscript = '/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/FACS_sorted_txMeth/R/CAMDAC_FACS_sorted_Samples.R'
ccdmp = load('/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/FACS_sorted_txMeth/data/CAMDAC_mt_vs_sorted_mt_tumour_vs_normal_DMP.RData')
ccdmp = copy(dt); rm(dt)
bdmp = load('/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/FACS_sorted_txMeth/data/normal_sorted_vs_bulk_DMP.RData')
bdmp = copy(tmp); rm(tmp)

# First, reproduce the 77% finding
## Data shows the percentage of FACS-mt DMPs (truth) that are also called by CAMDAC
ccfacs = ccdmp[ DMP_f_corr == T, ]
ccstats = ccfacs[, .(
    ncg = .N,
    ndmpb = sum(DMP_b_corr == T),
    ndmpt = sum(DMP_t_corr == T),
    fracb = sum(DMP_b_corr == T) / .N,
    fract = sum(DMP_t_corr == T) / .N
), by=.(patient, sample)]
aastats = ccfacs[, .(
    ncg = .N,
    ndmpb = sum(DMP_b_corr == T),
    ndmpt = sum(DMP_t_corr == T),
    fracb = sum(DMP_b_corr == T) / .N,
    fract = sum(DMP_t_corr == T) / .N
)]

ccstats
aastats

# Create a figure from ccstats that shows the spread for bulk and pure.
figdir = fs::dir_create('figures/240214_FACS')
fracstats = tidyr::pivot_longer(ccstats, cols=c('fracb', 'fract'), names_to='type', values_to='frac') %>%
    mutate(pid=paste0(patient, "_", sample))

pfrac = ggplot(fracstats, aes(x=pid, y=frac, fill=type)) +
    geom_hline(yintercept=mean(ccstats$fract), linetype='dashed', color='#ff93a5') +
    geom_hline(yintercept=mean(ccstats$fracb), linetype='dashed', color='#80daf8') +
    geom_bar(stat='identity', position='dodge') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title='FACS-mt DMPs called by CAMDAC', x='Patient', y='Fraction of DMPs') +
    scale_fill_manual(values=c('blue', 'red')) +
    geom_text(aes(label=ncg), y=0.9, size=3) +
    scale_y_continuous(n.breaks=20) 


ggsave(paste0(figdir, '/FACS-mt_DMPs_called_by_CAMDAC.png'), pfrac, width=10, height=5)


unique(ccdmp[, .(patient, sample)])

# All percentages point to >92%, even among bulk
sum(ccdmp$DMP_f == ccdmp$DMP_t) / nrow(ccdmp)
sum(ccdmp$DMP_f_corr == ccdmp$DMP_t_corr) / nrow(ccdmp)
sum(ccdmp$DMP_b == ccdmp$DMP_t) / nrow(ccdmp)
sum(ccdmp$DMP_b == ccdmp$DMP_f) / nrow(ccdmp)


ccdmp[1,]
