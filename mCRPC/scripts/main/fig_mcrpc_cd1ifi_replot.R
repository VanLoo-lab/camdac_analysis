# Load functions
source('src/mcrpc_rhmr.R')
options(bitmapType='cairo')

# Load data
dt = qs::qread('data/out/fig_mcrpc_cd1ifi.qs')

# Read objects from data
p_methl = dt[[1]]
p_dcdif = dt[[2]]
rm(dt); gc()

# Get QC pass samples
qcpass = fread('data/out/hmrs/meta_mcrpc.tsv') %>% 
    dplyr::filter(tPASS==TRUE) %>% dplyr::pull(sample_id)

# Load data
xdata = p_methl[[2]]$data %>% dplyr::filter(variable %in% qcpass)

# Generate the delta between bulk and pure
xlong = tidyr::pivot_wider(xdata, names_from='dataset', values_from='mean_meth',
    id_cols=c("variable", "cg_ix", "chrom", "start", "end")) %>% data.table() %>%
    dplyr::mutate(delta = pure-bulk)

# Cluster pure methylation rates and return sample identifiers
pmat = xlong  %>%
    dplyr::select(cg_ix, variable, delta) %>% 
    tidyr::pivot_wider(names_from=variable, values_from=delta) %>%
    dplyr::select( -cg_ix) %>% as.matrix()

# Cluster pmat by rows and then return the dendogram order of the column names
porder = hclust(dist(t(pmat)))$order
pnamesord = colnames(pmat)[porder]
pnamesord = pnamesord[pnamesord %in% qcpass]
dcord = dplyr::filter(xlong, variable %in% pnamesord) %>%
    dplyr::mutate(vfact = factor(variable, levels=pnamesord))

# Drop samples with >20% CpGs missing
cg_20pct = 0.1*max(dcord$cg_ix)
vcg = dcord %>% dplyr::group_by(vfact) %>% dplyr::summarise(nsum=sum(is.na(delta))) %>%
    dplyr::filter(nsum<8500) %>% pull(vfact)
dcord = dplyr::filter(dcord, vfact %in% vcg)

# Generate plot
pdelta = ggplot(dcord, aes(x=factor(cg_ix), y=vfact)) +
    geom_tile(aes(fill=delta)) +
    scale_fill_gradient2(low='blue', mid='white', high='red') +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank()) + 
    theme(axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
    labs(y="mCRPC tumours", x='CpGs', fill='Δ​​m') 
ggsave('figures/out/cd1ifi_delta_4by2.png', plot=pdelta, width=4, height =2)
ggsave('figures/out/cd1ifi_delta_4by2.pdf', plot=pdelta, width=4, height =2)
ggsave('figures/out/cd1ifi_delta_4by4.png', plot=pdelta, width=4, height =4)
ggsave('figures/out/cd1ifi_delta_4by4.pdf', plot=pdelta, width=4, height =4)
qs::qsave(pdelta, 'figures/out/p_cd1ifi_delta.qs')

# Get regions that end up in the plot and find where CD1 and IFI start
# As a function of the number of CpGs displayed.
#    For example, if CD1 stats CpG 50/100, this function prints (0.5).
#      Use this as a multiplier on the plot length to manually set start position
xreg = unique(dcord[, .(chrom, start, end, cg_ix)])

whereis <- function(xreg, chr, xstart, xend){
    # get the row indexes of the regions that contain the start and end
    # and return what fraction they are of the total rows
    # note: all chroms already in chr1
    xreg$ix = 1:nrow(xreg)
    minfrac = min(xreg[chrom==chr & start>=xstart]$cg_ix) / max(xreg$cg_ix) 
    maxfrac = max(xreg[chrom==chr & end<=xend]$cg_ix) / max(xreg$cg_ix) 
    print(paste0("CD1/IFI start: ", round(minfrac,3)))
    print(paste0("CD1/IFI end: ", round(maxfrac,3)))
}

# Where is cd1a
whereis(xreg, "1", 158254424, 158258269)

# Where is ifi16
whereis(xreg, "1", 159009907, 159055151)

# Print regions for anntation
xreg[cg_ix==min(xreg$cg_ix) | cg_ix==max(xreg$cg_ix)]
