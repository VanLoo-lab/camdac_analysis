# Load data and plot read length by phasing
library(ggplot2)
library(dplyr)
library(data.table)

phfiles = fs::dir_ls('data/source/PHASED', glob='*.qs')
phfiles
loadphased <- function(fid) {
    sid = stringr::str_extract(fs::path_file(fid), "LTX\\d+_R\\d+")
    ph = qs::qread(fid)
    # Select only CpGs above 10X
    ph = ph[nreads>=10]
    # Count the number of CpGs overlapping SNPs and overlapping SNVs by read length
    prl = ph[, .(nCpG = .N, oSNP=sum(npsnp>=3), oSNV=sum(npsnv>=3)), by=.(rlen)]
    # Calculate proportions for oSNP and oSNV
    prl[, poSNP:=oSNP/nCpG]
    prl[, poSNV:=oSNV/nCpG]
    prl$sample_id=sid
    prout = prl %>% tidyr::pivot_longer(c(poSNP, poSNV, oSNP, oSNV, nCpG), names_to="type", values_to="count")
    return(prout)
}

phased = data.table::rbindlist(lapply(phfiles, loadphased))
figdir = fs::dir_create('figures/240304_readlen_phased')
plt1 = ggplot(phased[type=="nCpG"], aes(x=rlen, y=count, color=sample_id)) + geom_line() + 
    geom_point() + theme_classic() + scale_color_brewer(palette="Set1")
ggsave(fs::path(figdir, "readlen_nCpG.png"), plt1, width=6, height=4)

plt2 = ggplot(phased[type=="oSNP"], aes(x=rlen, y=count, color=sample_id)) + #geom_line() + 
    geom_point() + theme_classic() + scale_color_brewer(palette="Set1") +
    labs(x="Read length", y="CpGs on reads overlapping SNPs")
ggsave(fs::path(figdir, "readlen_oSNP.png"), plt2, width=6, height=4)

plt3 = ggplot(phased[type=="poSNP"], aes(x=rlen, y=count, color=sample_id)) + #geom_line() + 
    geom_point() + theme_classic() + scale_color_brewer(palette="Set1") +
    labs(x="Read length", y="Proportion of CpGs (by read length)\n overlapping SNPs")
ggsave(fs::path(figdir, "readlen_poSNP.png"), plt3, width=6, height=4)

plt4 = ggplot(phased[type=="oSNV" & count>0], aes(x=rlen, y=count, color=sample_id)) + #geom_line() + 
    geom_jitter() + theme_classic() + scale_color_brewer(palette="Set1") +
    geom_smooth(method="lm", se=FALSE, color="black") + labs(x="Read length", y="CpGs on reads overlapping SNVs")
ggsave(fs::path(figdir, "readlen_oSNV.png"), plt4, width=6, height=4)

