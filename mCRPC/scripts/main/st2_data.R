# Aim: Collate purity and ploidy metrics for supplementary table 2

library(data.table)
library(dplyr)
library(stringr)
library(readr)

# Setup
rwes = fread('data/archive/2021_data/input_files/output_wes_overlap.txt')
rwgs = fread('data/archive/2021_data/input_files/210223_wgs_overlap_result.txt')
rmc = fread('data/pre/f2_mcrpc_cna/mcrpc_wgbs_wgs_overlap.csv')
rcrc = fread('data/source/CRC_WGBS_WGS_overlap.txt')
rcrr = fread('data/source/CRC_RRBS_WGS_overlap.txt')

# Request CRUKID fix
txmeta = fread('data/meta/tracerxPublicationKey_170221.txt')
txcruk = txmeta[[2]]; names(txcruk) = txmeta[[1]]

# Setup helper to convert CRUK IDs to tracerX
# Add cruk identifier
i_convert_tx_to_cruk = function(txid) {
    # Split input string
    txsplit = str_split(txid,'_')[[1]]

    crukid = txcruk[txsplit[1]]
    if(is.na(crukid)) {
        return(NA)
    } else {
        # Add number e.g. _R1
        return(paste0(crukid, '_', txsplit[2]))
    }
}
convert_tx_to_cruk = Vectorize(i_convert_tx_to_cruk)

# Collate

r0 = rwes[,
    .(
        cohort='NSCLC',
        txid=rrbs_sample,
    camdac_source='RRBS',
    dna_source='WES',
        camdac_purity=rrbs_purity,
        dna_purity=wes_purity,
        camdac_ploidy=rrbs_ploidy,
        dna_ploidy=wes_ploidy)
] %>% unique() %>% 
    dplyr::mutate(
        txid= txid %>% str_replace('\\.', '_') %>% str_replace('LTX', 'LTX0'),
    ) %>% dplyr::mutate(
        txid= convert_tx_to_cruk(txid)
    ) 


r1 = rwgs[,
    .(
    cohort='NSCLC',    
    txid=rrbs_sample,
    camdac_source='RRBS',
    dna_source='WGS',
        camdac_purity=rrbs_purity,
        dna_purity=wgs_purity,
        camdac_ploidy=rrbs_ploidy,
        dna_ploidy=wgs_ploidy)
] %>% unique() %>% 
    dplyr::mutate(
        txid= txid %>% str_replace('\\.', '_') %>% str_replace('LTX', 'LTX0'),
    ) %>% dplyr::mutate(
        txid= convert_tx_to_cruk(txid)
    ) 

r2 = rmc %>%
    dplyr::mutate(
        camdac_source='WGBS',
        dna_source='WGS',
        cohort='mCRPC'
    ) %>%
    dplyr::select(
        cohort,
        txid=sample_id,
        camdac_source,
        dna_source,
        camdac_purity=purity_b,
        dna_purity=purity_a,
        camdac_ploidy=ploidy_b,
        dna_ploidy=ploidy_a) 

r3 = rcrc %>%
    dplyr::mutate(
        camdac_source='WGBS',
        dna_source='WGS',
        cohort='CRC'
    ) %>%
    dplyr::select(
        cohort,
        txid=sample_id,
        camdac_source,
        dna_source,
        camdac_purity=purity_b,
        dna_purity=purity_a,
        camdac_ploidy=ploidy_b,
        dna_ploidy=ploidy_a)

r4 = rcrr %>%
    dplyr::mutate(
        camdac_source='RRBS',
        dna_source='WGS',
        cohort='CRC'
    ) %>%
    dplyr::select(
        cohort,
        txid=sample_id,
        camdac_source,
        dna_source,
        camdac_purity=purity_b,
        dna_purity=purity_a,
        camdac_ploidy=ploidy_b,
        dna_ploidy=ploidy_a)


# Combine all results
r = rbind(r0, r1, r2, r3, r4, fill=T)

data.table::fwrite(r, 'data/out/st2_data.csv')
