library(ggplot2)
library(stringr)
library(data.table)
library(patchwork)

get_sid <- function(x) str_extract(x, "DTB-\\d+")

# Functions ---
plt_snv_phased <- function(x, xname, yname, ptitle){
    ggplot(x, aes(x=get(xname), y=get(yname))) +
    # Add linear model line with no error bars
    geom_point(alpha=0.15, size=0.5) +
    geom_abline(intercept=0, slope=1, color="navyblue", linetype='dashed', alpha=0.5) +
    geom_smooth(method=lm, se=FALSE, color="red", alpha=0.5) +
    theme_classic() +
    scale_x_continuous(limits=c(0,1)) +
    scale_y_continuous(limits=c(0,1)) +
    labs(title=ptitle, x=xname, y=yname)
} 

# Load SNVmt data for each sample
load_snvmt_data <- function(curdir, snvdp, main) {
    sid = get_sid(curdir)
    svp = helper_load_snv_phased(curdir)
    svp$qname=NULL
    svp = unique(svp)
    
    # Load CAMDAC deconvolved data for the same sample
    pure_file = fs::dir_ls(fs::path(main, sid), recurse=T,glob="*pure.csv.gz")
    pure = fread(pure_file)
    pure = pure[, .(chrom, start, end, m, cov, mhdiw=m_x_high-m_x_low,
            m_t_raw, m_t, cov_t, mhdiw_t=m_t_high-m_t_low, nA, nB,
            purity_p=purity, ploidy_p=ploidy, seg_start, seg_end, CG_CN)]
    # Annotate pure with SNV mt
    setkey(pure, chrom, start, end)
    setkey(svp, chrom, start, end)
    svpure = merge(svp, pure, all.x=T, by=c("chrom", "start", "end"))
    
    # Load CAMDAC SNV VAF data
    snvdata = fread(fs::path(snvdp, paste0(sid, "_snv_dpc.txt")))
    setnames(snvdata, "pos", "hap_POS")
    setnames(snvdata, "ref", "snv_REF")
    setnames(snvdata, "alt", "snv_ALT")
    setnames(snvdata, "VAF", "snv_VAF")
    setnames(snvdata, "mut_cn", "snv_MUTCN")
    setnames(snvdata, "ccf", "snv_CCF")
    svpures = merge(svpure, snvdata, all.x=T, by=c("chrom", "hap_POS"))
    
    return(svpures)
}

# This is now a lot fo data, and only for one chromosome (19K CpGs).
#   Should we save intermediate data? (svpures)
#   For now, limit to all we would need for plotting
valphas_to_minimal <- function(svpures){
valphas = svpures[, .(chrom, start, end, hap_POS, hap_BAF, snv_REF, snv_ALT, snv_CCF, snv_MUTCN, snv_VAF,
    altm, atcm, refm, rtcm, major_cn, minor_cn, purity, ploidy, m, mhdiw, cov, m_t, mhdiw_t, cov_t, CG_CN)]
valphas = valphas[!is.na(altm) & !is.na(m_t)]
valphas = valphas[minor_cn==0]
return(valphas)
}

calc_absdiff <- function(x, xval, yval){
    xsplit = split(x, x$patient)

    xdiff = lapply(xsplit, function(x){
        metric_all = abs(x[[xval]] - x[[yval]])
        bstats = boxplot.stats(metric_all)$stats
        res = data.frame(t(bstats)) 
        names(res) = c("metric_min", "metric_q1", "metric_med", "metric_q3", "metric_max")
        res$metric_mean = mean(metric_all, na.rm=T)
        res$metric_sd = sd(metric_all, na.rm=T)
        res$metric_ins = paste0(xval, "_", yval)
        res$e_nrpcc = mean(x$atcm/(x$major_cn+x$minor_cn), na.rm=T)
        res$patient = x$patient[[1]]    
        return(res)
        })
    xd = rbindlist(xdiff)
    return(xd)
}

plt_ndiag <- function(x){
    pp1 = ggplot(x, aes(x=e_nrpcc, y=metric_med, color=metric_ins)) +
        geom_point() + 
        # Add a line with secodn order polynomial fit
        geom_smooth(method=lm, formula=y~poly(x,3), se=FALSE, alpha=0.5) +
        theme_classic()

    x$patient = forcats::fct_reorder(x$patient, x$nrpcc)
    pp2 = ggplot(x, aes(x=patient, y=metric_med, color=nrpcc)) +
        geom_point(aes(shape=metric_ins)) +
        # draw a line between the points of the groups defined by patient
        geom_line(aes(group=patient)) +
        scale_color_viridis_c() +
        theme_classic() +
        theme(axis.text.x = element_blank())

    pp_patch = patchwork::wrap_plots(pp1, pp2, nrow=1)
    return(pp_patch)
}
