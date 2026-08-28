library(data.table)
library(ggplot2)
library(tidyr)
options(bitmapType='cairo')

outdir = fs::dir_create("data/pre/gsea")
ora_res = data.table::fread(fs::path(outdir, 'ora_res.tsv'))
hm = dplyr::filter(ora_res, hallmark=="hall")

hmw = ora_res %>% 
    dplyr::mutate(lpadj2=-log10(padj2)) %>%
    dplyr::mutate(index=1:nrow(.)) %>%
    dplyr::select(pathway, dtype, hallmark, source, dset, lpadj2) %>%
    tidyr::pivot_wider(names_from=dset, values_from=lpadj2)
gp = ggplot(hmw, aes(x=bulk, y=pure, fill=hallmark)) + geom_point(shape=21, alpha=.6) + 
    facet_wrap(dtype~source, scales='free') + geom_abline(intercept=0, slope=1, linetype='dashed') + theme_minimal()
ggsave(fs::path(outdir, 'ora_res_hm.png'), gp, width=10, height=10, dpi=300)

