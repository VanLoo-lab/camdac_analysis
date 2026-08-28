
# Load libraries
source('src/mcrpc_rhmr.R')

figdir = fs::dir_create('figures/240417_rhmr_mcrpc')
hmrd = '/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/hmrs'
bhmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/bulk_hmr.qs'
phmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/pure_hmr.qs'
meta = fread('data/out/hmrs/meta_mcrpc.tsv')
#cgi = fread('data/meta/hg38_cgi_open_cln.csv.gz') 
qcpass = dplyr::filter(meta, tPASS==TRUE) %>% dplyr::pull(sample_id)
dscores = get_dscores() # Add immune labels
# /nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/mcrpc_rna/figures/danaher.png

# Load metho objects
metho = qs::qread("data/out/metho/cdif_metho.qs")
bmetho = flatten_metho(metho, 'bulk')
pmetho = flatten_metho(metho, 'pure')
methl = rbind(bmetho, pmetho)

# Annotate immune status
methl = merge(methl, dscores, by='variable', all.x=T)
methl = dplyr::filter(methl, variable %in% qcpass)
log_info("Loading plots")
#p_methl = plot_cdif(methl, qcpass)
log_info("Plotting data")
#ggsave(fs::path(figdir, 'cd1_ifi6_meth_cg.png'), p_methl, width=10, height=14, dpi=300)

# Get rhmr data
bhmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/bulk_hmr.qs'
phmr = '/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/hmrs/pure_hmr.qs'
brhmr = load_rhmrs(bhmr)
prhmr = load_rhmrs(phmr)


# 3) Show the HMRs
bcdif <- get_cdif_meth(brhmr) %>% dplyr::mutate(dataset='bulk')
pcdif <- get_cdif_meth(prhmr) %>% dplyr::mutate(dataset='pure')
dcdif <- rbind(bcdif, pcdif)
#p_dcdif <- plot_cdif(dcdif, qcpass)
#ggsave(fs::path(figdir, 'cd1_ifi6_meth.png'), p_dcdif, width=10, height=10, dpi=300)

# Save all important objects into a list for later
#qs::qsave(list(p_methl, p_dcdif), 'data/out/fig_mcrpc_cd1ifi.qs')

