# Plot intersection, purity and ploidy for solutions tables from CAMDAC
library(data.table)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
options(bitmapType='cairo')

# Functions -----

#' Return a plot dataframe for each of the three plots we wish to create:
#'  intersection | purity | ploidy
get_intersect_data <- function(seg_a, seg_b, name_a, name_b){
    seg_a = fread(seg_a)
    seg_b = fread(seg_b)
    # Validate inputs
    stopifnot(
        all(names(seg_a) == names(seg_b))
    )

    # Add names
    seg_a$name = name_a
    seg_b$name = name_b

    # Add suffix to relevant fields
    suffix_add = setdiff(names(seg_a), c("sample_id", "chrom"))
    for(i in suffix_add){
        setnames(seg_a, i, paste0(i, "_a"))
        setnames(seg_b, i, paste0(i, "_b"))
    }

    # Make CNA segments intersection plot data. Must have same sample IDs.
    sids = unique( intersect(seg_a$sample_id , seg_b$sample_id) )
    cna_i = lapply(sids, function(x) intersect_cna_segments(seg_a, seg_b, x))
    cna_c = data.table::rbindlist(lapply(cna_i, function(x) calculate_copy_overlap(x$combined)))
    
    return(
        list(
            res=cna_c,
            intersects=cna_i
        )
    )
}

plot_intersect_data <- function(data, outfile){
    p_purity = plot_purities(data)
    p_ploidy = plot_ploidies(data)
    p_intersect = plot_intersects(data)
    p_wrap = patchwork::wrap_plots(p_purity, p_ploidy, p_intersect, nrow=1)
    ggsave(outfile, p_wrap, width=12, height=4)
}

plot_intersects <- function(e){
    e = dplyr::arrange(e, desc(perc_overlap))
    e$samples = factor(seq(nrow(e)))
    x= tidyr::pivot_longer(e, c("perc_overlap", "perc_overlap_1"), names_to="overlap_type", values_to="perc_overlap")
    x$overlap_type = ifelse(x$overlap_type == "perc_overlap", "exact", "within_1")
    xax=paste0(x$name_a[1], "_", x$method_a[1])
    yax=paste0(x$name_b[1], "_", x$method_b[1])
    ggplot(x, aes(x=samples, y=perc_overlap, fill=overlap_type)) +
        geom_point(shape=23, color='black') +
        scale_fill_manual(values=c("darkgreen", "orange")) +
        labs(x="Sample", y="Percent overlap", title="Overlap") +
        theme_classic() +
        scale_y_continuous(limits=c(0,100)) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        # place legend within plot at bottom right
        theme(legend.position=c(0.25, 0.25))

}

plot_ploidies <- function(x){
    xax=paste0(x$name_a[1], "_", x$method_a[1])
    yax=paste0(x$name_b[1], "_", x$method_b[1])
    x$ploidy_a = ifelse(x$ploidy_a > 6, 6, x$ploidy_a)
    x$ploidy_b = ifelse(x$ploidy_b > 6, 6, x$ploidy_b)
    # Set last label on X and Y axis to '6+'
    xax_labels = seq(0,6,1)
    xax_labels[length(xax_labels)] = "6+"
    yax_labels = seq(0,6,1)
    yax_labels[length(yax_labels)] = "6+"
    ggplot(x, aes(x=ploidy_a, y=ploidy_b)) +
    geom_abline(intercept=0, slope=1, color='grey', linetype='dashed') +
    geom_abline(intercept=2, slope=1, color='grey', linetype='dashed') +
    geom_abline(intercept=-2, slope=1, color='grey', linetype='dashed') +
    geom_point(shape=23, color='black', fill='skyblue') +
    scale_x_continuous(limits=c(0,6), breaks=seq(0,6,1), labels=xax_labels) +
    scale_y_continuous(limits=c(0,6), breaks=seq(0,6,1), labels=yax_labels) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x=stringr::str_to_upper(xax), y=stringr::str_to_upper(yax), title="Ploidy") +
    theme_classic()
}

plot_purities <- function(x){
    xax=paste0(x$name_a[1], "_", x$method_a[1])
    yax=paste0(x$name_b[1], "_", x$method_b[1])
    ggplot(x, aes(x=purity_a, y=purity_b)) + 
    geom_abline(intercept=0, slope=1, color='grey', linetype='dashed') +
    geom_point(shape=23, color='black', fill='purple') +
    scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x=stringr::str_to_upper(xax), y=stringr::str_to_upper(yax), title="Purity") +
    theme_classic()
}

intersect_cna_segments <- function(ssa, ssb, sid){
    # ssx is a dataframe with the fields: chrom, start_x, end_x, major_cn_x, minor_cn_x, name_x, method_x
    # Warning: Must be applied to single samples only

    ssa = dplyr::filter(ssa, sample_id == sid)
    ssb = dplyr::filter(ssb, sample_id == sid)

    # Validate inputs
    is_missing = any(c(nrow(ssa) == 0, nrow(ssb) == 0))
    if(is_missing){
        return(list(combined=data.frame()))
    }

    # Create GRanges objects of both inputs
    gr_a = makeGRangesFromDataFrame(ssa, start.field="start_a", end.field="end_a", keep.extra.columns=TRUE)
    gr_b = makeGRangesFromDataFrame(ssb, start.field="start_b", end.field="end_b", keep.extra.columns=TRUE)

    # Get the intersection between the two frames and keep mcols
    intersection = pintersect(findOverlapPairs(gr_a, gr_b))
    meta_a <- mergeByOverlaps(gr_a, intersection)
    int_a <- meta_a$intersection
    mcols(int_a) <- mcols(meta_a$gr_a)

    # >>>
    meta_b <- mergeByOverlaps(gr_b, intersection)
    int_b <- meta_b$intersection
    mcols(int_b) <- mcols(meta_b$gr_b)

    # >>>
    int_a <- data.frame(int_a)
    int_b <- data.frame(int_b)
    combined <- merge(int_a, int_b, by=c("seqnames","start","end"))

    return(
        list(
            combined=combined,
            meta_a=meta_a,
            meta_b=meta_b
        )
    )

}

calculate_copy_overlap <- function(combined){
    # Validate
    if(nrow(combined)==0){
        return(data.frame())
    }
  # Return the percentage of bases at segment intersect with same CN state across assays
  copy_overlap_bool <- (combined$major_cn_a == combined$major_cn_b) & (combined$minor_cn_a == combined$minor_cn_b)
  major_bool <- (combined$major_cn_a == combined$major_cn_b)
  # Limit to regions of overlap
  bases_copy_overlap <- sum(combined[copy_overlap_bool,]$width.x)
  bases_total_intersection <- sum(combined$width.x)
  perc_copy_overlap <- bases_copy_overlap * 100 / bases_total_intersection

  # As above but for within_1
    major_within_1 <- (abs(combined$major_cn_a - combined$major_cn_b) <= 1) & (combined$minor_cn_a == combined$minor_cn_b)
    minor_within_1 <- (combined$major_cn_a == combined$major_cn_b) & (abs(combined$minor_cn_a - combined$minor_cn_b) <= 1)
    copy_overlap_bool_1 <- major_within_1 | minor_within_1
    bases_copy_overlap_1 <- sum(combined[copy_overlap_bool_1,]$width.x)
    bases_total_intersection_1 <- sum(combined$width.x)
    perc_copy_overlap_1 <- bases_copy_overlap_1 * 100 / bases_total_intersection_1
  
  res = data.frame(
    sample_id=combined$sample_id.x[1],
    purity_a=combined$purity_a[1],
    purity_b=combined$purity_b[1],
    ploidy_a=combined$ploidy_a[1],
    ploidy_b=combined$ploidy_b[1],
    method_a=combined$method_a[1],
    method_b=combined$method_b[1],
    name_a=combined$name_a[1],
    name_b=combined$name_b[1],
    bases_overlap=bases_copy_overlap[1],
    bases_intersection=bases_total_intersection[1],
    perc_overlap=perc_copy_overlap[1],
    bases_overlap_1=bases_copy_overlap_1[1],
    bases_intersection_1=bases_total_intersection_1[1],
    perc_overlap_1=perc_copy_overlap_1[1]
  )
  return(res)
}


# Main -----

f_sega = data.frame(
    sample_id="S1",
    purity=.5,
    ploidy=2.2,
    goodness=80,
    dataset="wgs",
    method="ascat",
    chrom=1,
    start=sample(1:1000, 20)
) %>% mutate(end=start+200)
data.table::fwrite(f_sega, "test.csv")

f_segb = f_sega[1:5,] %>%
    mutate(
        dataset="wgbs",
        method="ascat"
    )
data.table::fwrite(f_segb, "test2.csv")

# Matcherd ASCAT solutions
data = get_intersect_data("test.csv", "test2.csv", "wgs", "wgbs")
plot_intersect_data(data$res, "test.png")
