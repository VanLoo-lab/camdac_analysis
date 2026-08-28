# Unclear why, but WGBS DMR annotatios do not contain the promoter regions as expected.
## To fix, generate again from ENSMBL
## Code adapted from original at : /nemo/project/proj-vanloo-secure/working/nmensah/projects/wgbs_camdac_dev/analysis/archive/create_camdac_refs/create_camdac_references.Rmd
library(data.table)
library(biomaRt)
library(GenomicRanges)
library(IRanges)
library(dplyr)
outdir = fs::dir_create('data/archive/2024_data/wgbs_prom_regions')
mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#mart_f =  useEnsembl(biomart="ENSEMBL_MART_FUNCGEN", dataset="hsapiens_regulatory_feature") 


gene_names <- getBM(
  attributes=c(
    "chromosome_name","start_position","end_position","strand",
    "external_gene_name","ensembl_gene_id",
    "transcription_start_site","transcript_start","transcript_end",
    "ensembl_transcript_id","transcript_gencode_basic"
  ), 
    filters = "chromosome_name", values = c(1:22, "X", "Y"), mart = mart
  )
ogene_names=data.table::copy(gene_names)
gene_names = data.table(gene_names)
setnames(gene_names, "chromosome_name", "chrom")
setnames(gene_names, "start_position", "start")
setnames(gene_names, "end_position", "end")
qs::qsave(gene_names, fs::path(outdir, 'gene_names.qs'))

# Promoter annotation ----
# Limit to gencode transcripts as Lili did
gencode = gene_names[!is.na(transcript_gencode_basic) & transcript_gencode_basic!=""]
gencode$strand = fifelse(gencode$strand==1, "+", "-")
TSSs <- GenomicRanges::GRanges(
  seqnames=gencode$chrom,
  strand=gencode$strand,
  ranges=IRanges(start=gencode$transcription_start_site, 
                 end = gencode$transcription_start_site)
)
# Lili's annotations used 2.5kb and 250b. I'm using GRanges defaults which are more conservative at 
# 2kb and 200b. We create the promoter gene names in the same way
promoters <- GenomicRanges::promoters(TSSs, upstream=2000, downstream=200)
promoters <- data.table(data.frame(promoters)) %>% dplyr::select(chrom=seqnames, start, end)
promoters$external_gene_name = gencode$external_gene_name
promoters$ensembl_gene_id = gencode$ensembl_gene_id
promoters[, gene_prom_ids:=paste(
  external_gene_name, ensembl_gene_id, sep=","
)]
promoters[, gene_prom_ids:=stringr::str_replace(gene_prom_ids, "^,", "")] # Clean starts

# Overlap promoters but collapse names by unique identifier
setkey(promoters, chrom,start,end)
qs::qsave(promoters, fs::path(outdir, 'promoters.qs'))


# Add overlap from existing annotation
