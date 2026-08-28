# From the CAMDAC pipeline files, calculate the percentage of 1000G SNPs covered by CpG alleles in our reference dataset
setwd('~/projects/paper_camdac')
library(data.table)
library(dplyr)
library(GenomicRanges)

# Define loader function
load_plfs = function(path){
    plfs = fs::dir_ls(path)
    loci = GRanges()
    for (x in plfs){
        d_ = load(x)
        loci = c(loci, loci_subset)
    }
    loci = data.table(data.frame(loci))
    return(loci)
}

# For RRBS
plfs = load_plfs('~/pvs/pipeline_files/rrbs/pipeline_files/loci_files/')
r_issnp = !is.na(plfs$POS)
r_iscg = ifelse(plfs$width > 1, T, F)
sum(r_issnp & r_iscg) / sum(r_issnp) # 49.2% of SNP alleles covered by CpG alles
sum(r_issnp & r_iscg) / sum(r_iscg) # 25.6% of CpG alleles covered by SNP alleles

# For WGBS
wlfs = load_plfs('~/pvs/pipeline_files/wgbs/hg38/loci_files/')
wlfs = wlfs[!duplicated(wlfs[,.(seqnames, start, end)]),]
## Select ramwas CpG set
## http://www.people.vcu.edu/~ashabalin/RaMWAS/cpgset_ref_hg38_chr1-22.rds
cgset = readRDS('data/meta/cpgset_ref_hg38_chr1-22.rds')
cgset_df = do.call(rbind, lapply(names(cgset), function(chr) {
    data.table(seqnames = chr, start = cgset[[chr]], end = cgset[[chr]]+1) # Assuming end is start + 1
}))
cgset_df$iscg = T
wlfs = merge(wlfs, cgset_df, by=c('seqnames', 'start', 'end'), all.x=T)
w_issnp = !is.na(wlfs$POS) # 27m SNP
w_iscg = !is.na(wlfs$iscg) & wlfs$iscg
sum(w_issnp & w_iscg) / sum(w_issnp) # 16.3% of SNP alleles covered by CpG alles
sum(w_issnp & w_iscg) / sum(w_iscg) # 16.9% of CpG alleles covered by SNP alleles

# Add 1000G
g1000 = fs::dir_ls('~/proj-vanloo/reference_files/human/references/1000genomes/1000genomes_2014_phase3_loci/', glob='*genomesloci*_chr*.txt')
gdf = lapply(g1000, function(x){
    dt = fread(x); names(dt)= c("seqnames", "pos"); return(dt)
})
gdf = data.table::rbindlist(gdf) # 84m positions
gdf$seqnames = paste0('chr', gdf$seqnames)
gdf$olstart=T
gdf$start=gdf$pos
cgnew = merge(cgset_df, gdf, by=c('seqnames', 'start'), all.x=T)
gdf$olend=T;gdf$olstart=NULL;gdf$end=gdf$start;gdf$start=NULL;
cgnew = merge(cgnew, gdf, by=c('seqnames','end'), all.x=T)
cgnew$is_g1000 = (!is.na(cgnew$olstart) & cgnew$olstart) | (!is.na(cgnew$olend)  & cgnew$olend)
sum(cgnew$is_g1000)*100/nrow(cgnew) # 5.86% of CpG alleles covered by 1000G snps; Ignore and stick with reference files

# However, the question is wrong. We want the % of CGs covered by 1000G SNPs.