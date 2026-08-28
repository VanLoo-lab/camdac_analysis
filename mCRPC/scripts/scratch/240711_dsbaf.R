# Plot the WGS downsampled BAF
# Original: /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/scripts/archive/231002_cpu_wgs_down_generate.R\

# Setup
library(qs)
library(ggplot2)
library(dplyr)
library(data.table)
options(bitmapType='cairo')

# Load data
tnl = qs::qread('/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/archive/2023_data/230916_wgbs_ascat/wgs_ascat_ds_inputs/cache/DTB-104/DTB-104_tnl.qs')

# Adapted function to plot just the three BAF values

flip_n_hets <- function(vec, n_to_change){
    # Get indices of FALSE values
    false_indices <- which(!vec)

    # Randomly sample from false_indices
    indices_to_change <- sample(false_indices, n_to_change)

    # Set the chosen indices to TRUE
    vec[indices_to_change] <- TRUE
    return(vec)
}


# probabilistic approach to assign heterozygous SNPs directly from tumour BAF profiles
is_het <- function(x, y, pbin = 0.01, probHom = .99, na.rm = TRUE) {
  flag <- logical(length = length(x))
  flag <- !(pbinom(unlist(x), size = unlist(y), prob = probHom, log.p = FALSE) > pbin)
  flag <- ifelse(flag == TRUE, "Heterozygous", "Homozygous")
  return(flag)
}


ascat_wgs_ds_baf <- function(tnl){

    tnl$ix = seq(1, nrow(tnl))
    
    # Assign tumor-only hets
    baf_set <- ifelse(tnl$new_baf < .5, (1-tnl$new_baf)* tnl$new_depth, tnl$new_baf * tnl$new_depth)
    hets <- is_het(baf_set, tnl$new_depth) == "Heterozygous"
    hets <- flip_n_hets(hets, 20)
    to_het = tnl[ hets ]

    # Assign tumor-normal hets
    tn_het = tnl[ new_nBAF > 0.2 & new_nBAF < 0.8  ]
    to_het$ix = seq(1, nrow(to_het))

    # Set WGBS BAF
    pwgbs = ggplot(to_het, aes(x=ix, y=BAF)) + geom_point(alpha=0.01, size=0.05, color='black') + theme_bw() + labs(x="SNP-index", title="DTB-104 : WGBS", y="BAF")+
  theme(plot.title = element_text(hjust = 0.5))
    # Set original WGS BAF
    porig = ggplot(to_het, aes(x=ix, y=wgbaf)) + geom_point(alpha=0.01, size=0.05, color='red') + theme_bw() + labs(x="SNP-index", title="DTB-104 : WGS", y="BAF")+
  theme(plot.title = element_text(hjust = 0.5))
    # Set WGS Downsampled BAF
    pwgds = ggplot(to_het, aes(x=ix, y=new_baf)) + geom_point(alpha=0.01, size=0.05, color='red') + theme_bw() + labs(x="SNP-index", title="DTB-104 : Downsampled WGS", y="BAF") +
  theme(plot.title = element_text(hjust = 0.5))

    pwrap = patchwork::wrap_plots(
        list(porig, pwgbs, pwgds),
        ncol=1, nrow=3
    )
    return(pwrap)
}

bafds <- ascat_wgs_ds_baf(tnl)
figdir = fs::dir_create("/nemo/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/figures/240711_exitsem_104ds")
figout = fs::path(figdir, 'wgs_ds_baf.png')
ggsave(figout, bafds, width=12, height=7, dpi=200)


