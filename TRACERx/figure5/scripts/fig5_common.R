# Shared helpers for the Figure 5 SNV-phased DMR panels (5A pies, 5B-D examples).
# Classification and colours are those of Notes_phase_methylation_SNV_v2.R (scripts/, L885-1146).

## colours: Notes L1126-1131 (black/grey40/olivedrab/royalblue4/hotpink4 + lighter hypo) as they
## appear in the published panels (fills were drawn with alpha; sampled here as solid colours)
snv_dmr_cols <- c(both_hyper = "#101010", both_hypo = "#9D9D9D", cis_hyper = "#7A9342", cis_hypo = "#C5DF8E",
                  LOH_hyper = "#394B8C", LOH_hypo = "#8B9FE7", trans_hyper = "#884A69", trans_hypo = "#F4A5CD",
                  none = "#FEFEFA")

## Classify SNV-phased bins relative to the SNV (mutant) allele (Notes L885-1021).
## dt: table from SNV_DMRs_v2.RData. Returns the phased bins (>= 3 phased CpGs) with DMR_AS/DMR_1/DMR_2 and `flag`.
classify_snv_dmrs <- function(dt, n_AS_DMPs = 3L) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  suppressWarnings(dt[, c("m_n_DMPs", "m_t_DMPs") := lapply(.SD, as.numeric), .SDcols = c("m_n_DMPs", "m_t_DMPs")])
  dt[, DMR_AS := ifelse(!is.na(consec_AS_DMPs) & consec_AS_DMPs >= n_AS_DMPs, "DMR", as.character(NA))]
  dt[, DMR_1 := ifelse(!is.na(consec_DMPs_1) & consec_DMPs_1 >= n_AS_DMPs, "DMR", as.character(NA))]
  dt[, DMR_2 := ifelse(!is.na(consec_DMPs_2) & consec_DMPs_2 >= n_AS_DMPs, "DMR", as.character(NA))]
  dt <- dt[CpG_phased_counts >= 3, ]                                     # phased bins (>= 3 phased CpGs)
  
  dt[, flag := ifelse(!is.na(DMR)&CN_1<=0.5&!is.na(m_t_2_DMPs)&m_t_2_DMPs-m_n_DMPs>0, "LOH_hyper",
               ifelse(!is.na(DMR)&CN_1<=0.5&!is.na(m_t_2_DMPs)&m_t_2_DMPs-m_n_DMPs<0, "LOH_hypo",
               ifelse(!is.na(DMR)&is.na(DMR_AS)&CN_1>0.5&
                      !is.na(m_t_1_DMPs)&m_t_1_DMPs-m_n_DMPs>(0.2)&
                      !is.na(m_t_2_DMPs)&m_t_2_DMPs-m_n_DMPs>(0.2), "both_hyper",
               ifelse(!is.na(DMR)&is.na(DMR_AS)&CN_1>0.5&
                      !is.na(m_t_1_DMPs)&m_t_1_DMPs-m_n_DMPs<(-0.2)&
                      !is.na(m_t_2_DMPs)&m_t_2_DMPs-m_n_DMPs<(-0.2), "both_hypo",
               ifelse(CN_1>0.5&!is.na(DMR)&!is.na(DMR_AS)&!is.na(m_t_2_AS_DMPs)&
                      abs(m_t_2_AS_DMPs-m_n_DMPs)>abs(m_t_1_AS_DMPs-m_n_DMPs)&
                      m_t_2_AS_DMPs-m_n_DMPs>0, "cis_hyper",
               ifelse(CN_1>0.5&!is.na(DMR)&!is.na(DMR_AS)&!is.na(m_t_2_AS_DMPs)&
                      abs(m_t_2_AS_DMPs-m_n_DMPs)>abs(m_t_1_AS_DMPs-m_n_DMPs)&
                      m_t_2_AS_DMPs-m_n_DMPs<0, "cis_hypo",
               ifelse(CN_1>0.5&!is.na(DMR)&!is.na(DMR_AS)&!is.na(m_t_2_AS_DMPs)&
                      abs(m_t_2_AS_DMPs-m_n_DMPs)<abs(m_t_1_AS_DMPs-m_n_DMPs)&
                      m_t_1_AS_DMPs-m_n_DMPs>0, "trans_hyper",
               ifelse(CN_1>0.5&!is.na(DMR)&!is.na(DMR_AS)&!is.na(m_t_2_AS_DMPs)&
                      abs(m_t_2_AS_DMPs-m_n_DMPs)<abs(m_t_1_AS_DMPs-m_n_DMPs)&
                      m_t_1_AS_DMPs-m_n_DMPs<0, "trans_hypo",
               "none"))))))))]
  dt
}
