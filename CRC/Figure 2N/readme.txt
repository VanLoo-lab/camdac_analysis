WGBS_new: use phased reads to get VAF (so as long as the SNV has a CpG phased to it on the same read, it's included.) BAF1 is extracted from VCF files directly; BAf2 is calculated from WGS REF/ALT counts. Directly generated from T.asm_hap_stats.csv.gz.
WGBS_latest: reformatted to match RRBS objects; no filtering on depth is applied. 



CpGs phased to 88832/2983 SNVs from WGBS (value1/value2 stands for all sites / clonal LOH and filtered by MCN sites). See script in WGBS_new: Figure_2N_WGBS.RDS is generated with *asm_hap_stats.csv.gz - Summary statistics for each phased SNV.
Notice that the count is SNV, not CpG. Multiple CpGs could be phased to the same SNV.