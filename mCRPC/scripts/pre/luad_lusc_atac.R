# LUAD and LUSC - ATAC-seq genes
# Source: https://elifesciences.org/articles/62522/figures#content
# Supplementary Figure 2 : https://iiif.elifesciences.org/lax/62522%2Felife-62522-fig1-figsupp2-v2.tif/full/,1500/0/default.jpg
# Selected genes distinguishing subtypes
source('src/at2_basal_regions.R')

# Inputs ----
rb_root = '/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/txMeth_DMPs'
outdir = fs::dir_create('data/pre/at2_basal')

# Load dataset with AT2 and Basal gene regions
data <- load_wang_data()
areg <- ranges_and_proms(data)

# Load annotation for samples
sanno <- load_sample_anno()

# Load DNA methylation per sample
rb_files = fs::dir_ls(rb_root, recurse=T, glob='*purified_tumour_DMP_calls.RData')

# Loop over methylation files and run processing function
doParallel::registerDoParallel(cores=10)
res = foreach(ff=rb_files) %dopar% {
    
    # Get specific sample data
    sdata <- parse_sdata(ff, sanno)
    if(nrow(sdata)==0){return(NULL)}

    # Load df_mt
    load(ff); x=df_mt
    if(!("DMP_t_corr" %in% names(x))){return(NULL)}

    # Parse over aregions
    oldata = process_meth_dmp(x, areg)
    resdata = cbind(oldata, sdata)
    return(resdata)
}
doParallel::stopImplicitCluster()
qs::qsave(res, fs::path(outdir, 'at2_basal_regions.qs'))
res = data.table::rbindlist(res)
qs::qsave(res, fs::path(outdir, 'at2_basal_regions.qs'))