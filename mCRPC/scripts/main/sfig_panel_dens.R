# Aim: Plot panel of normal performance for each tissue type, showing the density, DMP percentages and ratio of hyper and hypomethylation

options(bitmapType='cairo')
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
library(logger)
source('src/sfig_panel.R')
figdir = fs::dir_create("figures/240305_panels")

## The boxplot data are in these objects under object [[1]]
## The densities are in these objects, under object [[2]]
valfiles = fs::dir_ls("./data/out/panel/validation/")

# x = './data/out/panel/validation/test_object_male.qs'
# test = list(x,x,x,x,x,x,x)
val = purrr::discard(valfiles, ~stringr::str_detect(.x, 'test'))

# Validation plot -----
## /camp/project/proj-vanloo-secure/working/nmensah/projects/paper_camdac/data/out/panel/validation/wgbs_bone_normal_male.qs

cplots = lapply(val, function(e){
  ptitle = stringr::str_split(fs::path_file(e), "_", simplify=T)[,2]
  psex = stringr::str_split(fs::path_file(e), "_", simplify=T)[,4]
  ptitle = paste0(ptitle, ":", psex)
  ptitle = stringr::str_to_upper(fs::path_ext_remove(ptitle))
  pipe_cplots(e, ptitle) 
})
cfinal = patchwork::wrap_plots(cplots, nrow=1)
qs::qsave(cfinal, '~/projects/paper_camdac/data/out/panel/cfinal.qs')
ggsave(fs::path(figdir, "panel_validation.png"), cfinal, width=10, height=5)
ggsave(fs::path(figdir, "panel_validation.pdf"), cfinal, width=10, height=5)

log_info("Panel Validation Complete")

# # Saturation plot -----

satdir = "./data/out/panel/saturation/"
satfiles = fs::dir_ls(satdir)
satd = rbindlist(lapply(satfiles, process_sat))
satd = dplyr::mutate(satd,
  sex = stringr::str_split(fname, "_", simplify=T)[,4],
  tissue = stringr::str_split(fname, "_", simplify=T)[,2]
  )
satd = dplyr::filter(satd, tissue!='object')
satmed = data.table(satd)[, .(mn=max(current_n)), by='fname'][mn>=10,]$fname
satmed = dplyr::filter(satd, fname %in% satmed)
plt = plot_sat_med_line(satmed)
figfile = fs::path(figdir, "normal_panel_sat.png")
figfile_pdf = fs::path(figdir, "normal_panel_sat.pdf")
ggsave(figfile, plt, width=5, height=3)
ggsave(figfile_pdf, plt, width=5, height=3)

qs::qsave(
  list(satd, cfinal, plt),
  'data/out/panel/figdata.qs'
)

# New validation plot : Figure S6D ----

# Create a figure showing the boxplot across all val files
# First, combine datasets
parse_box_data <- function(x){
  dt = qs::qread(x)
  vald = process_val(dt[[1]]) %>% dplyr::mutate(
    source=fs::path_ext_remove(fs::path_file(x))
  )
  return(vald)
}

loaded = parallel::mclapply(val, parse_box_data)
loaded = data.table::rbindlist(loaded)

# Clean source
loaded$source = loaded$source %>%
  stringr::str_replace_all('wgbs_', '') %>%
  stringr::str_replace_all('normal_', '') %>%
  stringr::str_replace_all('_', '\n')

lclean = loaded %>% dplyr::filter(DMP!='perc_ns')
lns = loaded %>% dplyr::filter(DMP=='perc_ns')

# Write functino to generate boxplot
make_dmp_all <- function(vald){
  vald$DMPtype = ifelse(vald$DMP=='perc_ns', 'N.S.', 'DMP')
  ddmp = unique(vald[ , .(n_cg, perc=sum(percentage)), by=c('item', 'source', 'DMPtype')])
  dboxplot = ggplot(ddmp, aes(x=source, y=perc)) +
    theme_bw() +
    geom_boxplot(aes(color=DMPtype), outlier.shape=NA) +
    scale_fill_manual(values=c('red', 'grey', 'black')) +
    scale_color_manual(values=c('red', 'grey', 'black')) +
    scale_y_continuous(breaks=seq(0,30,5)) +
    geom_jitter(aes(group=DMPtype, fill=DMPtype), alpha=0.5, size=0.5, color='black', position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.75)) +
    labs(x='Tissue', y='CpGs (%)') +
    facet_wrap(~DMPtype, ncol=1, scales='free_y') +
    theme(
      strip.background = element_blank(),
      legend.position = 'none'
    ) +
    scale_y_continuous(n.breaks=5)
    #+
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

dplts = make_dmp_all(loaded)
ggsave(fs::path(figdir, 'val_boxplot_all.pdf'), dplts, width=8, height=3.4)
data.table::fwrite(loaded, fs::path(figdir, 'val_boxplot_all.csv'))


# Count DMP information

ddmp = unique(lclean[ , .(n_cg, perc=sum(percentage)), by=c('item', 'source')])
mean(ddmp$perc)
sd(ddmp$perc)
# > mean(ddmp$perc)
# [1] 0.6296341
# > sd(ddmp$perc)
# [1] 0.8611067

dlns = unique(loaded[ , .(n_cg, perc=sum(percentage)), by=c('item', 'source')])
mean(dlns$perc)
sd(dlns$perc)
# > mean(dlns$perc)
# [1] 8.943415
# > sd(dlns$perc)
# [1] 6.759003