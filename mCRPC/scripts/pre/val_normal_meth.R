#!/camp/apps/eb/software/R/3.6.2-foss-2019b/bin/Rscript
#SBATCH --job-name=CAMDAC_MAIN
#SBATCH --output=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/files/logs/val_%A_%a.out
#SBATCH --error=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/files/logs/val_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=300G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Create a file listing locations of normal methylation and assocaited metadata
library(data.table)
suppressPackageStartupMessages(library(dplyr))
source("src/camdac_metho.R")

## Setup outputs
outdir = fs::dir_create('data/pre/val_normal_meth')

## Normal methylation for mCRPC samples
wnman = fread("data/source/manifest_wgbs_ac.tsv")
wnman = wnman[ type=="normal" ]

## SKIP: Add normal methylation from Chinese atlas
mch = fread("data/source/manifest_wgbs_ac_li.tsv")
mch = mch[ type=="normal" & tissue=="prostate" ]
mch$tissue="prostate_ch"
wnman = rbind(wnman, mch)

## TODO: Run with tracerX normal allele counts files
mtx = fread("data/source/manifest_wgbs_ac_tracerx_normal.tsv")
wnman = rbind(wnman, mtx)

## Write manifest to file
fwrite(wnman, "data/source/manifest_wgbs_ac_normal.tsv", sep="\t")
stopifnot(all(sapply(wnman$ac_file, fs::file_exists)))

## Split methylation data for tissues with at least 3 samples
wnman = wnman[tissue %in% (dplyr::count(wnman, tissue, sex_anno) %>% dplyr::filter(n>=3) %>% pull(tissue) %>% unique()),]
wnsplit = split(wnman, list(wnman$tissue, wnman$sex_anno))
wnsplit = wnsplit[sapply(wnsplit, function(x) nrow(x)>0)] 

## Create a methylation matrix object for each dataset table
## Requires fields : dataset, sample_id, ac_file, tissue, sex_anno
meth_names = names(wnsplit)
meth_names = purrr::discard(meth_names, str_detect(meth_names, 'prostate_ch'))

for(mn in meth_names){
    manifest = wnsplit[[mn]]
    outfile = manifest[1,] %>% dplyr::transmute(
        outname=paste0("wgbs_", tissue, "_", type, "_", sex_anno, ".qs")
    ) %>% pull(outname)
    outpath = fs::path(outdir, outfile)

    if(fs::file_exists(outpath)){
        log_info("Skipping {mn}")
        next
    }
    log_info("Processing {mn}")
    meth_object = load_camdac_metho(manifest, chunksize=20, ncores=10)

    # Fix pos of non-CAMDAC counts (chr prostate. Doesn't work for others)
    if(manifest$tissue[[1]] == "prostate_ch"){
        meth_object$pos$start = meth_object$pos$end -1
    }

    qs::qsave(meth_object, outpath)
    log_info("Saved {mn}")
}
