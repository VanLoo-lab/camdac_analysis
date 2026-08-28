# Cell scoring method from Danaher 2017
## The gene symbols from the original manuscript did not work,
## so instead I am leveraging gene identifiers from Lili's personal dataset.
## Here, we recreate the clusters on the WGBS data

# # Provides GEX score of TIL for each sample
# # https://jitc.bmj.com/content/5/1/18
# # https://static-content.springer.com/esm/art%3A10.1186%2Fs40425-017-0215-8/MediaObjects/40425_2017_215_MOESM2_ESM.pdf

# Data transferred from nemo
setwd('/nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/mcrpc_rna')
renv::load("~/projects/wgbs_camdac_dev")
library(data.table)
library(dplyr)
library(ggplot2)
library(janitor)
library(org.Hs.eg.db)

# Load RNA data
rna = fread('data/mcrpc_rna_cpm_symbol.tsv')
genes = fread('data/combined_immune_gene_set_with_entrezids.txt', header=T)
sum(genes$gene %in% rna$gene)/nrow(genes) # 97% overlap

# Log-transform RNA CPM matrix
log_expr = function(x, min=0.5){ifelse(x > min, log(x), log(min))}
rcols = setdiff(names(rna), 'gene')
grna = rna[, ..rcols]
trna <- grna[, lapply(.SD, log_expr), .SDcols = rcols]
trna$gene = rna$gene

# First, set all the rows to have the same genes in the grna matrix
stopifnot(nrow(rna) == nrow(grna))
gsel = which(rna$gene %in% genes$gene & !duplicated(rna$gene))
filtrna = rna[gsel,]
filtgrna = as.matrix(grna[gsel,])
filtgenes = filtrna$gene


make_matrix <- function(genes, filtgenes, fname="cell_type"){
    genes$i = 1
    fields = c("gene", fname, "i")
    x = tidyr::pivot_wider(genes[, ..fields], id=fname, names_from="gene", values_from="i",
        values_fn=function(x) ifelse(any(x)==1, 1, 0))
    x[is.na(x)] = 0
    x = as.matrix(x)
    rnamesx = x[,1]
    x = x[,-1]
    x = apply(x, 2, as.numeric)
    rownames(x) = rnamesx
    x = x[, filtgenes]
    return(x)
}

lin_mat = make_matrix(genes, filtgenes, "lineage")
cell_mat = make_matrix(genes, filtgenes, "cell_type")
cell_mat[1,]

dim(filtgrna)
dim(cell_mat)
dim(lin_mat)

lin_av = lin_mat %*% filtgrna
cell_av = cell_mat %*% filtgrna
stopifnot(all(colnames(lin_av) == colnames(cell_av)))
ests = rbind(lin_av, cell_av)

qs::qsave(ests, 'data/mcrcp_danaher_scores_ests.qs')
ests = qs::qread('data/mcrcp_danaher_scores_ests.qs')

make_hm <- function(ests, scale=T){
    if(scale){
        #ests = scale(ests)
        ests = t(scale(t(ests)))
    }
    ddist <- dist(t(ests))
    hc <- hclust(ddist, method = "ward.D")
    # Cut the tree in a sensible manner
    hc_clusters <- factor(cutree(tree = hc, k = 4))
    hc_anno = data.frame(
        clusters=hc_clusters,
        samples=names(hc_clusters)
    )
    stopifnot(all(hc_anno$samples == colnames(ests)))

    # Plot a heatmap of ests
    library(circlize)
    library(ComplexHeatmap)
    dev.off()
    ha <- HeatmapAnnotation(df=hc_anno, show_legend=F)
    phm = Heatmap(ests, color = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation=ha)
    res = list('heatmap'=phm, 'clusters'=hc_clusters, 'samples'=hc_anno$samples, 'data'=ests)
    return(res)
}

# Save to png in figures/240219_danaher_heatmap.png
options(bitmapType='cairo')

# draw heatmap
library(png)
library(gridExtra)
phm = make_hm(ests, T)
hmap = phm$heatmap
gb = grid.grabExpr(draw(phm$heatmap))
ggsave('figures/danaher_v2.png', gb, width=10, height=10, units="in")

qs::qsave(list(
    "rna"=ests,
    "scaled"=phm$data,
    "clusters"=phm$clusters,
    "heatmap"=phm$heatmap
), 'data/mcrcp_danaher_scores.qs')
