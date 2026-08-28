library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

motif1 <- "CCGG"
motif2 <- "CG"
hg38 <- BSgenome.Hsapiens.UCSC.hg38

matches <- lapply(seqnames(hg38)[1:24], function(seqname) {
    print(seqname)
    seq <- getSeq(hg38, seqname)
    match_ccgg <- matchPattern(motif1, seq)
    match_cg <- matchPattern(motif2, seq)
    
    ccgg_df <- data.frame(seqname = seqname, start = start(match_ccgg), end = end(match_ccgg), motif = motif1)
    cg_df <- data.frame(seqname = seqname, start = start(match_cg), end = end(match_cg), motif = motif2)
    
    rbind(ccgg_df, cg_df)
})

# Combine the results into a single data frame
matches_df <- do.call(rbind, matches)
names(matches_df)[[1]] = 'chrom'
data.table::fwrite(matches_df, 'data/meta/hg38_ccgg_cg_positions.csv.gz')
