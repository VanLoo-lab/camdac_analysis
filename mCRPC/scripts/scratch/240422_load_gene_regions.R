# Get annotations for gggenes for regions of interest
#  molecule gene start   end  strand orientation
#1  Genome1 genA 15389 17299 reverse           2

# Region for CD1IFI6
cchrom=1
cstart=156834461
cend=160189117

library("biomaRt")
listMarts(host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)


regions<-paste0(cchrom,":",cstart,":",cend)
attributes<-c("hgnc_symbol", "chromosome_name","start_position", "end_position", "strand", "band")
results=getBM(attributes=attributes,filters="chromosomal_region",values=regions,mart=ensembl)

