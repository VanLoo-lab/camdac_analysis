#!/camp/apps/eb/software/R/3.6.2-foss-2019b/bin/Rscript
#SBATCH --job-name=CAMDAC_MAIN
#SBATCH --output=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/files/logs/main_%A_%a.out
#SBATCH --error=/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/files/logs/main_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Create a matrix with the methylation values from all camdac pure files
source("src/camdac_metho.R")
library(data.table)
library(foreach)
library(doParallel)
library(fs)
library(logger)
library(qs)
library(stringr)

# For each file, make a metho object with bulk and pure data
pure_files = fs::dir_ls('data/source/main/', glob="*pure.csv.gz", recurse=T)

# Read all CAMDAC data as meth objectcs
doParallel::registerDoParallel(12)
metho_l = foreach(i=pure_files, .combine='c') %dopar% {
    pure_file = i
    sample_id = stringr::str_extract(pure_file, "DTB-[0-9]+")
    tmp = tempfile(fileext=".qs")
    log_info("Processing {sample_id}")
    metho = dispatch_read_metho(sample_id, pure_file)
    qs::qsave(metho, tmp)
    return(tmp)
}
doParallel::stopImplicitCluster()

# Combine objects into a single one
log_info("Combining {length(metho_l)} metho objects...")
metho = qs::qread(metho_l[[1]])
remainder = metho_l[-1]
for(i in seq_along(remainder)){
    log_info("Merging {i+1}...")
    obj = qs::qread(remainder[[i]])
    metho = merge_metho(metho, obj)
}

# Saving result
outdir = fs::dir_create("data/out/metho")
qs::qsave(metho, fs::path(outdir, "main_pure.qs"))
log_info("Done")

