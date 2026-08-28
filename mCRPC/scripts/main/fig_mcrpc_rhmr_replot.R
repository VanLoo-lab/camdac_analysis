# Aim: Plot rHMR heatmaps and CGI counts between bulk and pure tumours

# Previous work: 
#   sfig: /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/scripts/main/sfig_rHMRs.R
#   prepro: /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/scripts/pre/format_hmrs_matrix.R

# Load libraries
source('src/mcrpc_rhmr.R')
options(bitmapType='cairo')

# Inputs ----
figdir = fs::dir_create('figures/240417_rhmr_mcrpc_v2')
#hmrd = '/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/hmrs'
bhmr = 'data/out/hmrs/bulk_hmr.qs'
phmr = 'data/out/hmrs/pure_hmr.qs'
meta = fread('data/out/hmrs/meta_mcrpc.tsv')
cgi = fread('data/meta/hg38_cgi_open_cln.csv.gz') 
qcpass = dplyr::filter(meta, tPASS==TRUE) %>% dplyr::pull(sample_id)

# Get rhmr data
brhmr = load_rhmrs(bhmr)
prhmr = load_rhmrs(phmr)

# Select top N rHMRs
brhtop = get_top_rhmrs(brhmr, 5000)
prhtop = get_top_rhmrs(prhmr, 5000)

# Load CpG and heatmap annotations (CGI status)
brcg = build_cg_map(bhmr, meta, cgi)
prcg = build_cg_map(phmr, meta, cgi)

# Annotate pure regions with those that are in bulk and vice versa
hmr_sets <- find_rhmr_sets(brhmr, prhmr, brhtop, prhtop)
bhmr_sets = hmr_sets[dataset=='bulk']
phmr_sets = hmr_sets[dataset=='pure']

# Plot pure and bulk heatmaps
set.seed(76) # For different colours
bulk_figout = fs::path(figdir, 'bulk_rhmr_heatmap.png')
bulk_hm <- plot_rhmr_heatmap(brhtop, bhmr_sets, meta, brcg$hmr_cgi, bulk_figout)
pure_figout = fs::path(figdir, 'pure_rhmr_heatmap.png')
pure_hm <- plot_rhmr_heatmap(prhtop, phmr_sets, meta, prcg$hmr_cgi, pure_figout)

# Plot CGI counts datajj
# 1) More HMRs in pure than in bulk
hcounts = rbind(
    dplyr::mutate(brcg$bdata, dataset='bulk'),
    dplyr::mutate(prcg$bdata, dataset='pure')
) %>% dplyr::mutate(dataset = stringr::str_to_title(dataset))
qs::qsave(hcounts, fs::path(figdir, 'mcrpc_rhmr.qs'))
p_hcount = plot_hcounts(brcg, prcg)
ggsave(fs::path(figdir, 'hmr_counts_4by2.pdf'), plot=p_hcount, width=4, height=2)
ggsave(fs::path(figdir, 'hmr_counts_4by4.pdf'), p_hcount, width=4, height=4)
ggsave(fs::path(figdir, 'hmr_counts_2by2.pdf'), p_hcount, width=2, height=2)

# 2) Greater proportion of rHMRs in open seas in pure
bhm_scgi = count_bp_cgi(brhmr, cgi) %>% dplyr::mutate(dataset='bulk')
phm_scgi = count_bp_cgi(prhmr, cgi) %>% dplyr::mutate(dataset='pure')
scgi = rbind(bhm_scgi, phm_scgi) %>% dplyr::mutate(dataset=stringr::str_to_title(dataset))
p_scgi <- plot_cgi_counts(scgi)
# Save in a few dimension ranges
ggsave(fs::path(figdir, 'cgi_counts_4by2.pdf'), p_scgi, width=4, height=2)
ggsave(fs::path(figdir, 'cgi_counts_4by4.pdf'), p_scgi, width=4, height=4)
ggsave(fs::path(figdir, 'cgi_counts_2by2.pdf'), p_scgi, width=2, height=2)


# Generate CGI plot
p_cgplot = patchwork::wrap_plots(list(p_hcount, p_scgi), ncol=1)
ggsave(fs::path(figdir, 'cgi_counts.png'), p_cgplot, width=10, height=10, dpi=300)

# Statistical test: Wilcoxon signed-rank test for pairs
hlong = hcounts %>% 
    dplyr::select(sample_id, dataset, hmrs) %>%
    tidyr::pivot_wider(names_from=dataset, values_from=hmrs) %>%
    dplyr::filter(complete.cases(.))
wx = hlong$Bulk; wy = hlong$Pure
pval = wilcox.test(wy, wx, conf.int=T, paired=TRUE) # 1.04e-14
pval
# Caluclate median difference
diff = median(wy - wx); diff

# Statistical test: proportion of HMRs in open seas
sclong = scgi %>%
    dplyr::select(samp, dataset, prop_osea) %>%
    tidyr::pivot_wider(names_from=dataset, values_from=prop_osea) %>%
    dplyr::filter(complete.cases(.))
wx = sclong$Bulk; wy = sclong$Pure
pval = wilcox.test(wy, wx, conf.int=T, paired=TRUE) # 1.04e-14
pval
# Caluclate median difference
diff = median(wy - wx); diff

