library(data.table)

# Get CAMDAC SNP type from SNP positions in the reference dataset
camloci = fs::dir_ls('~/pvs/pipeline_files/wgbs/hg38/loci_files/')
xdata = data.table::rbindlist(
    lapply(camloci, function(x){
        print(x)
        lv = load(x)
        seqlevelsStyle(loci_subset) = 'NCBI'
        dt = as.data.frame(loci_subset) %>% data.table
        dtp = dt[!is.na(POS)]
        dtp[, CPOS:=paste0(seqnames, ':', POS)]
        dtp[, ref_alt:=paste0(ref, '>', alt)]
        dtp = dtp[!duplicated(dtp$CPOS)]
        dtp$Chromosome = dtp$seqnames
        dtp$Position = dtp$POS
        anno = dtp[, .(Chromosome, Position, ref_alt)]
        return(anno)
    })
)

# Reclassift ref_alt to one of the four complementary strands
xdata[, ref_alt_comp:=fcase(
    ref_alt == 'G>T', 'C>A',
    ref_alt == 'G>C', 'C>G',
    ref_alt == 'G>A', 'C>T',
    ref_alt == 'A>T', 'T>A',
    ref_alt == 'T>C', 'C>T',
    ref_alt == 'A>G', 'C>T',
    ref_alt == 'A>C', 'C>A',
    ref_alt == 'T>G', 'C>A',
    rep(TRUE, length(ref_alt)), ref_alt
)]

data.table::fwrite(xdata, 'data/meta/hg38_camdac_snp_ref_alt.csv.gz')

