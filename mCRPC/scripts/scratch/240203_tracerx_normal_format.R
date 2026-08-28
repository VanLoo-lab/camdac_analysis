# Format tracerX normals to CAMDAC as manifest for use in the panel validation pipeline

library(data.table)
library(readr)
library(dplyr)
mtx = readr::read_lines("data/source/tx_normal_meths.txt")
tsid = stringr::str_extract(mtx, "LTX\\d+")
tsid = paste0(tsid, "_N")
tdataset = "tracerx"

txmeta = fread('/camp/project/proj-swanton-vanloo/working/tx421/dataset/RRBS_bismark/post_inventory_qc_ascat_22_05_19.tsv') %>% 
    dplyr::select(tracerx_id, region_id, sex, tx421, PASS) %>% 
    dplyr::filter(is.na(PASS) & region_id=="N" & tx421==F) %>%
    dplyr::mutate(sid=stringr::str_replace(tracerx_id, "LTX0", "LTX")) %>%
    dplyr::mutate(sid = paste0(sid, "_N")) %>%
    dplyr::select(sample_id=sid, sex_anno=sex) %>%
    dplyr::mutate(sex_anno=stringr::str_to_lower(sex_anno))
txmeta

txdata = data.table(
    dataset="su2c",
    sample_id=tsid,
    ac_file=mtx,
    tissue="lung",
    type="normal"
)
txdata = merge(txdata, txmeta, by="sample_id")

data.table::fwrite(txdata, "data/source/manifest_wgbs_ac_tracerx_normal.tsv", sep="\t")
