
# prostate immune assignments are here:
## > TODO. Calculate from CIBERSORT absolute values: /nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/mcrpc_rna/scripts/assign_immune_status.R 

# prostate SNV_wt results are here:
#   phasing in these folders: /nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/results/snvph/DTB-071/SNV_phase/T
#   function to load is: helper_load_snv_phased(phase_dir, chrom="all", cln_output=T)
#   code is: /nemo/project/proj-vanloo-secure/working/nmensah/data/240104_SU2C_CAMDAC/code/snv_mt_pipeline_from_asm.R
#   here, SNV_wt is refm. We can filter for known tumor areas using the hap_CCF and minor_cn==0

# TODO:
# Assign immune status using Danaher gene list from lili.
## Running in mcrpc dataset with run_danaher_2017.R