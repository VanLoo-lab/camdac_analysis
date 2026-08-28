#!/camp/apps/eb/software/R/3.6.2-foss-2019b/bin/Rscript
#SBATCH --job-name=wgbs_asm_dtb
#SBATCH --output=/nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/wgbs_asm/DTB_077_asm_%j.out
#SBATCH --error=/nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/wgbs_asm/DTB_077_asm_%j.out
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Dependencies ----
# Setup for run on HMEM or CPU
# renv::load("~/secure/working/nmensah/projects/wgbs_camdac_dev")
options(warn = -1)
renv::load("/nemo/project/proj-vanloo/secure/working/nmensah/cpu_renv/cpu_camdac/")
Sys.setenv("RENV_PATHS_CACHE"="/nemo/project/proj-vanloo/secure/working/nmensah/cpu_renv")
suppressPackageStartupMessages(devtools::load_all("~/secure/working/nmensah/projects/VanLoo_CAMDAC/"))

# Call CAMDAC for DTB-003
config <- CamConfig(
    refs = "/nemo/project/proj-vanloo/secure/working/nmensah/pipeline_files",
    outdir = "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_SU2C_CAMDAC/results/240418",
    bsseq = "wgbs",
    lib = "pe",
    build = "hg38",
    n_cores = 10
)

tumor <- CamSample(id = "T", sex = "XY", patient_id = "DTB-077", bam = "/nemo/project/proj-vanloo/secure/working/nmensah/su2c_ac/bam_symlink/DTB-077_T.WGBS.bam")
germline <- CamSample(patient_id = "DTB-077", id = "N", sex = "XY", bam = "/nemo/project/proj-vanloo-secure/working/nmensah/data/210828_SU2C_normals/bam_normal/DTB-077-Pro-NT-BSlib-U12_S40_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.sorted.bam")
infiltrates <- CamSample(patient_id = "DTB-077", id = "N_I", sex = "XY", bam = "/camp/home/mensahn/data/210505_su2c_bams/DTB-077-Pro-NT-BSlib-U12_S40_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.sorted.bam")
origin <- CamSample(patient_id = "DTB-077", id = "N_O", sex = "XY", bam = "/camp/home/mensahn/data/220819_PON_ANALYSIS/panels/bam/SRR6156035_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam")

# Add ASM CNA caller
attach_output(tumor, config, "asm_cna", "/nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/main/DTB-077/Copynumber/T/DTB-077.T.cna.txt")
#attach_output(tumor, config, "counts", "/nemo/project/proj-vanloo/secure/working/nmensah/su2c_ac/camdac_symlink/DTB-077_T.SNPs.CpGs.all.sorted.csv.gz")
#attach_output(germline, config, "counts", "/nemo/project/proj-vanloo-secure/working/nmensah/data/230919_SU2C_ASCAT_MATCHED/results/DTB-077/Allelecounts/N/DTB-077.N.SNPs.CpGs.all.sorted.csv.gz")
attach_output(infiltrates, config, "counts", "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male/bone.m.csv.gz")
attach_output(origin, config, "meth", "/nemo/project/proj-vanloo-secure/working/nmensah/data/230403_WGBS_PANELS/male/ch_prostate_norm.m.csv.gz")

asm_pipeline(tumor, germline, infiltrates, origin, config)
#asm_hap=qs::qread("/camp/home/mensahn/asm_hap.qs");asm_meth_cna=qs::qread("/camp/home/mensahn/asm_meth_cna.qs");asm_meth=qs::qread("/camp/home/mensahn/asm_meth.qs");phase_map=qs::qread("/camp/home/mensahn/phase_map.qs")
