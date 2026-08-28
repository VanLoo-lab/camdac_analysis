# Aim: Investigate the effect of read length on CpG-SNP phasing in CAMDAC

# Background:
#   - NGMerge was run on LTX064 samples
#       : /nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/ltx_nugen_ngmerge.sh
#   - NGMerge BAM output was processed with CAMDAC allele counter taking SNV and SNP loci
#       : /camp/home/mensahn/data/240212_TRACERX_RAW/ltx_nugen_count_rrbs_reads.R (see sourced file for main read phasing function)
#   - The output of the allele counter is saved in the PHASED/*.qs files
#       : /nemo/project/proj-swanton-vanloo/working/TRACERx_ONT/UCL_SERVER_TX/ltx_nugen/PHASED

# Setup
library(ggplot2)
library(dplyr)
library(data.table)
library(qs)
library(patchwork)
library(GenomicRanges)
options(bitmapType='cairo')

# Read phased data results files. Each record is a CpG:read-length pair with summary statistics:
#   - rlen: read length
#   - nreads: Number of reads at CpG with given rlen
#   - npsnp: Number of reads overlapping CpG and a SNP site
#   - npsnv: Number of reads overlapping CpG and an SNV site
#   - nreads_non1k: Number of reads !is.na(POS) CpGs
#   - nreads_non1ksnp: Number of non1k reads overlapping SNPs
#   - nreads_non1ksnv: Number of non1k reads overlapping SNVs
phfiles = fs::dir_ls('data/source/PHASED', glob='*.qs')

# Load phased data
x = qread(phfiles[1]) # LTX064_R1

# Aggregate counts at CpG-overlapping reads
x = x[, .(records=.N, sum_reads=sum(nreads), sum_snp=sum(npsnp), sum_snv=sum(npsnv)), by=.(rlen)]
x[, percSNP:=(sum_snp*100)/sum_reads]

# Vis. peaks are at 67, 133, 201
dplyr::arrange(x[rlen>60], -sum_reads)[1:10,]
dplyr::arrange(x[rlen>140], -sum_reads)[1:10,]

# Show reads overlapping CpGs at each read length
p_o1 = ggplot(x, aes(x=rlen, y=sum_reads)) + geom_col() + theme_classic() + labs(y='Read counts', x='Read length')

# Show percentage of (potential) SNP-phased reads overlapping CpGs
p_02 = ggplot(x, aes(x=rlen, y=percSNP)) + geom_col() + theme_classic() + labs(y='Reads overlapping common\n 1000 Genomes SNP loci (%)', x='Read length')

# Join patched
plt = patchwork::wrap_plots(p_o1, p_02, nrow=1)
figdir = fs::dir_create('figures/241126_ngmerge_phased')
ggsave(fs::path(figdir, "readlen_phased.png"), plt, width=10, height=4)

# Sambamba view commands to observe IGV peaks
# sambamba view -f bam -F "sequence_length > 66 and sequence_length < 69" 064_R1_100MB.sorted.bam  > 66_filtered.bam
# sambamba view -f bam -F "sequence_length > 130 and sequence_length < 135" 064_R1_100MB.sorted.bam  > 133_filtered.bam
# sambamba view -f bam -F "sequence_length > 190 and sequence_length < 204" 064_R1_100MB.sorted.bam  > 201_filtered.bam

# Get BED file of het SNPs from referenc ein the 100MB region
lf = load('/camp/home/mensahn/pvs/pipeline_files/rrbs/pipeline_files/loci_files/loci.hg19.1.RData')
lf = loci_subset[start(loci_subset) < 100000000]
lf = lf[!is.na(lf$POS)]
rbed = data.table(
    chrom=as.character(seqnames(lf)),
    start=lf$POS,
    end=lf$POS+1
)
data.table::fwrite(rbed, fs::path(figdir, 'hetSNPs.bed'), sep='\t', col.names=FALSE)