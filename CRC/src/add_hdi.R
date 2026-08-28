sel_asm_hdi_pass <- function (x, allele, hdi_fields) {
    col_prefix <- ifelse(allele == "ref", "ref_", "alt_")
    fields <- paste0(col_prefix, hdi_fields)
    complete.cases(x[, ..fields])
}


calculate_asm_hdi_bulk <- function(meth_c, allele, itersplit=1e5){
    col_prefix <- ifelse(allele == "ref", "ref_", "alt_")
    counts <- paste0(col_prefix, "total_counts_m")
    meth <- paste0(col_prefix, "m")
        # Generate inputs for HDI calculation
    M <- round(meth_c[[counts]] * meth_c[[meth]], 0)
    UM <- meth_c[[counts]] - M
    # Calculate ASM hdi for bulk data
    u_hdi <- CAMDAC:::unique_calculate_counts_hdi(M, UM, n_cores = n_cores, itersplit = itersplit)
    u_hdi <- round(u_hdi, digits = 5)
    names(u_hdi) = c("M", "UM", paste0(col_prefix, "m_lo"), paste0(col_prefix, "m_hi"))
    # Map back to original M and UM values
    inpdf = data.frame(M=M, UM=UM)
    u_inpdf = merge(inpdf, u_hdi, by=c("M", "UM"), all.x=TRUE, sort=FALSE)
    return(u_inpdf[, 3:4])
}
assignInNamespace("calculate_asm_hdi_bulk", calculate_asm_hdi_bulk, ns="CAMDAC")
