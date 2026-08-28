
parse_sim_hdi_long <- function(x, i){
    # x: table with RRBS simulation results
    # i: index of simulation for reporting

    # Extract information for T1
    dt_1 = data.table(
        DMP_index = 1,
        DMP_true = x$DMP_true_status_t1n,
        DMP_AS_DM = x$DMP_t1n_AS_DM,
        mt_hdi_low = x$low_mt1,
        mt_hdi_high = x$high_mt1,
        mt_hdi_width = pmin(abs(x$high_mt1 - x$low_mt1), 1),
        mt = x$m_t1,
        mt_cov = x$cov_t1,
        mb_hdi_low = x$low_mb1,
        mb_hdi_high = x$high_mb1,
        mb_hdi_width = pmin(abs(x$high_mb1 - x$low_mb1), 1),
        mb = x$m_b1,
        mb_cov = x$cov_b1,
        mn = x$m_n1,
        mn_cov = x$cov_n1,
        purity = x$p1,
        CN = x$CN_1
    )

    # Extract information for CG_2
    dt_2 = data.table(
        DMP_index = 2,
        DMP_true = x$DMP_true_status_t2n,
        DMP_AS_DM = x$DMP_t2n_AS_DM,
        mt_hdi_low = x$low_mt2,
        mt_hdi_high = x$high_mt2,
        mt_hdi_width = pmin(abs(x$high_mt2 - x$low_mt2), 1),
        mt = x$m_t2,
        mt_cov = x$cov_t2,
        mb_hdi_low = x$low_mb2,
        mb_hdi_high = x$high_mb2,
        mb_hdi_width = pmin(abs(x$high_mb2 - x$low_mb2), 1),
        mb = x$m_b2,
        mb_cov = x$cov_b2,
        mn = x$m_n2,
        mn_cov = x$cov_n2,
        purity = x$p2,
        CN = x$CN_2
    )

    # Combine extracted data
    dt = rbind(dt_1, dt_2)
    dt = cbind(data.table(sim = i), dt)

    # Return data.table
    return(dt)
}