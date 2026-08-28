# Aim: Repeat FACS validation from the manuscript, stratifying by immune cell type.
## This should be complemented by SNV_wt methylation
renv::load("~/projects/wgbs_camdac_dev")
library(data.table)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(fs)
library(qs)
library(foreach)
library(doParallel)
library(logger)
source('src/snvmt_functions.R')
source('src/val_snvmt.R')
devtools::load_all("~/projects/VanLoo_CAMDAC")
options(bitmapType='cairo')
get_sid <- function(x) str_extract(x, "DTB-\\d+")

qcfile = fs::path("data/out/qc/mcrpc_tumor.tsv")
outdir = dir_create("data/main/f3_immune/mcrpc")
outdir_root = fs::path_dir(outdir)
figdir = dir_create("figures/220221_f3_immune")
snvph= "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/snvph"
snvdp = "/nemo/project/proj-vanloo/secure/working/nmensah/wgs_cna_run/results/snv_dpc_txt" 
main= "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/main"

# FACS data are here:
# /nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/FACS_sorted_txMeth/data

# Phased SNVs are in here:
#txph = fs::dir_ls("/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX", recurse=T, glob="*SNV_phased*.RData")
#ex = load(txph[1])
#df_phased_methylation[1,]
# N/B: Lili has offered to do this for the RRBS. Leave and move on to DMR analysis.

# Tx immune assignments are taken from RDS where Lili had Rosenthal paper data
#imm = fread('/nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/immune/tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt')

# Load danaher scores and assignments
da = qs::qread('/nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/mcrpc_rna/data/mcrcp_danaher_scores.qs')

# Generate table of immune assignments
immasn = data.table(clust=da$clusters, sample_id=names(da$clusters))
coldid = da$clusters['DTB-024']
hotid = da$clusters['DTB-190']
immasn[, `:=`(immune=dplyr::case_when(
    clust==coldid ~ "cold",
    clust==hotid ~ "hot",
))]
immasn = immasn[!is.na(immune)]

# Load the SNV_wt data for all samples

# Main ----

# Get files with SNV phasing data
inval = fs::dir_ls('/camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/pre/f4_snvmt/results/', glob="*.csv.gz")
# Filter by QC pass samples
qcres = fread(qcfile)[tPASS==T]$patient
inval = inval[match(qcres,str_extract(inval, "DTB-\\d+"))]

# Get a vector with all the tissue panels
meta = fread('data/out/s1/table_s1.csv')
panels = list(
    bone=fread('/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male/bone.m.csv.gz'),
    liver=fread('/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male/liver.m.csv.gz'),
    lymph=fread('/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male/lymph.m.csv.gz')
)

# Load snvmt data for current sample
doParallel::registerDoParallel(cores=10)
res = foreach(curval=inval, .errorhandling="pass", .combine='rbind') %dopar% {

    log_info("Processing {curval}")
    sid = get_sid(curval)
    valpha = fread(curval)
    valpha$sample_id = sid

    # Get the correspondin tissue panel for this sample
    panid = meta[patient_id==sid & sample_id=="T",]$sample_site
    tpanel = panels[[panid]]
    stpanel = tpanel[, .(chrom, start, end, m_n=m, cov_n=cov, m_n_low=m_x_low, m_n_high=m_x_high, panel=panid)]

    # Overlap
    valpha$chrom = as.character(valpha$chrom)
    stpanel$chrom = as.character(stpanel$chrom)
    setkey(valpha, chrom, start, end)
    setkey(stpanel, chrom, start, end)
    valpan = merge(valpha, stpanel, all.x=T, by=c("chrom", "start", "end"))
    
    # Add immune status
    valpan = merge(valpan, immasn, all.x=T, by="sample_id")
    
    # Call DMPs
    mbdiff = valpan$m_n - valpan$refm
    M_n = round(valpan$m_n * valpan$cov_n,0)
    UM_n = valpan$cov_n - M_n
    M = round(valpan$refm * valpan$rtcm,0)
    UM = valpan$rtcm - M
    pdmp <- dmp_call_pipe(mbdiff, M_n, UM_n, M, UM,
        mtdiff = NULL,
        effect_size = 0.2, prob = 0.99, itersplit = 1e6, ncores=1
    )
    valdm = cbind(valpan, pdmp)
    return(valdm)
}

# Save object
outobj = fs::path(outdir, "immune_dmp.qs")
qs::qsave(res, outobj)

# Generate plot comparing immune hot vs cold mdiff distributions
# Note: DMPs can't be used
res = qs::qread(outobj)

library(ggplot2)
res = res[!is.na(res$mbdiff) & !is.na(res$immune),]
op = ggplot(res, aes(x=mbdiff, fill=immune)) + geom_density(alpha=0.5) + theme_minimal() + labs(title="Immune DMPs", x="Methylation difference", y="Density") + theme(legend.position="top") +
    facet_wrap(~panel, nrow=1) + theme(legend.position="top") +scale_fill_manual(values=c("hot"="red", "cold"="blue"))
op
# t test res on the mbdiff for hot vs cold
tt = t.test(res$mbdiff[res$immune=="hot"], res$mbdiff[res$immune=="cold"])
# get pvalue 
pval = round(tt$p.value,3)
options(bitmapType="cairo")
figdir = fs::dir_create("figures/240221_f3_immune")
ggsave(fs::path(figdir, "mcrcp_immune_snvwt.png"), op, width=12,  height=6) +
    labs(title=paste0("MCRCP SNV_wt vs Panel (p=", pval, ")"))

# Same again but do histogram
op2 = ggplot(res, aes(x=mbdiff, fill=immune)) + geom_histogram(alpha=0.5, position="identity", bins=50) + theme_minimal() + labs(title="Immune DMPs", x="Methylation difference", y="Density") + theme(legend.position="top") +
    facet_wrap(~panel, nrow=1) + theme(legend.position="top") +scale_fill_manual(values=c("hot"="red", "cold"="blue")) +
    scale_y_continuous(n.breaks=20)
op2
ggsave(fs::path(figdir, "mcrcp_immune_snvwt_hist.png"), op2, width=12,  height=6) +
    labs(title=paste0("MCRCP SNV_wt vs Panel (p=", pval, ")"))
# Filter to those that passed QC

# Show the methylation residuals plot
library(dplyr)
res[, .(mean_diff=mean(mbdiff), phased_cg=.N),by=c('panel', 'immune')] %>% dplyr::arrange(panel, immune)
