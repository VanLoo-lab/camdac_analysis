# Aim: Plot rHMR heatmaps and CGI counts between bulk and pure tumours

# Previous work: 
#   sfig: /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/scripts/main/sfig_rHMRs.R
#   prepro: /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/scripts/pre/format_hmrs_matrix.R

# Load libraries
setwd("~/projects/paper_camdac")
source('src/mcrpc_rhmr.R')

# Inputs ----
figdir = fs::dir_create('figures/241008_rhmr_stats') 
hmrd = '/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/hmrs'
bhmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/bulk_hmr.qs'
phmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/pure_hmr.qs'
meta = fread('data/out/hmrs/meta_mcrpc.tsv')
cgi = fread('data/meta/hg38_cgi_open_cln.csv.gz') 
qcpass = dplyr::filter(meta, tPASS==TRUE) %>% dplyr::pull(sample_id)

# Plot CGI counts datajj
# 1) More HMRs in pure than in bulk
hcounts = qs::qread('data/out/hmrs/hmr_counts.qs')

p_hcount = ggplot(hcounts, aes(x=dataset, y=hmrs, color=dataset, fill=dataset)) + 
    geom_jitter(width=0.2, shape=21, color='black') + 
    geom_violin(alpha=0.7) + 
    theme_bw() +
    labs(color="", fill="", x="", y="Hypomethylated regions from MethlySeekR\n(100bp windows)") + 
    theme_bw() + 
    scale_y_log10(n.breaks=10, labels=label_number(scale_cut=cut_si(unit=""))) + 
    scale_fill_manual(values=c('lightgrey', '#5aa3f4')) +
    scale_color_manual(values=c('lightgrey', '#5aa3f4')) +
    theme(legend.position='none')


# 2) Greater proportion of rHMRs in open seas in pure
scgi = qs::qread('data/out/hmrs/scgi_counts.qs')
p_scgi <- ggplot(scgi, aes(x=dataset, y=prop_osea, fill=dataset)) +
    geom_jitter(width=0.2, shape=21, color='black') + 
    geom_boxplot(alpha=0.7) +
    theme_minimal() + scale_y_continuous(limits=c(0,1)) + scale_fill_manual(values=c('lightgrey', '#5aa3f4')) + 
    scale_color_manual(values=c('lightgrey', '#5aa3f4')) +
    theme(legend.position='none') + labs(y="Proportion of HMRs in Open Seas", x="") + theme_bw() + theme(legend.position="none")

# Generate CGI plot
p_cgplot = patchwork::wrap_plots(list(p_hcount, p_scgi), ncol=1) 
ggsave(fs::path(figdir, 'cgi_counts.png'), p_cgplot, width=3, height=4, dpi=300)
ggsave(fs::path(figdir, 'cgi_counts.pdf'), p_cgplot, width=3, height=4, dpi=300)

qs::qsave(list(p_hcount, p_scgi), fs::path(figdir, 'cgi_counts.qs'))
