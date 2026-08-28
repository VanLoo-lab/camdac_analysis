
process_val <- function(val) {
  val_long <- val %>%
    tidyr::pivot_longer(cols =  c("perc_hyper", "perc_hypo", "perc_ns"), names_to = "DMP", values_to = "percentage") %>%
    data.table()
  val_long$cg_label = paste0(round(val_long$n_cg/1e6,2), "m")
  return(val_long)
}

plot_val_long <- function(val_long) {
  ggplot(val_long, aes(y=item, x=percentage, fill=DMP)) + geom_col(color='black') +
    theme_bw() + scale_fill_manual(
      labels=c("perc_hyper"="Hyper", "perc_hypo"="Hypo", "perc_ns"="N.S."),
      values=c("Pink", "Skyblue", "white")) +
    labs(y="", x="Percentage CpG sites") +
    geom_text(aes(label=cg_label, x=max(percentage)+0.5, y=item), size=3) +
     facet_grid(rows=vars(tissue, sex), scales="free_y", space="free_y") +
  theme(strip.background = element_rect(fill="white"), strip.text = element_text(face="bold"))

}

process_sat <- function(satfile) {
  sat = qs::qread(satfile)
  if(nrow(sat) == 0) return(NULL)
  satmed = sat %>%
    group_by(item, current_n) %>%
    summarise(med=median(100-perc_non), iqr=IQR(100-perc_non)) %>%
    data.table()
  satmed$fname = fs::path_file(satfile) %>% fs::path_ext_remove()
  return(satmed)
}

plot_sat_med <- function(satmed){
    ggplot(satmed, aes(x=factor(current_n), y=med, color=fname)) +
    geom_line(aes(group=item), position=position_dodge(.5)) +
    geom_linerange(aes(ymin=med-iqr, ymax=med+iqr),  position=position_dodge(.5)) +
    theme_bw() + geom_point(position=position_dodge(.5)) +
    labs(x="Samples in panel", y= "Median DMPs across comparisons (%) ") +
    theme(legend.position='bottom') 
}

plot_sat_med_line <- function(satmed){
    # Count patients and adjust filenames accordingly
    ncomp = satmed[, .(item, fname)] %>% unique() %>% dplyr::count(fname) 
    ncomp$fname_new = paste0(ncomp$fname, " (n=", ncomp$n, ")")
    satmed = merge(satmed, ncomp, by='fname', all.x=T)
    satmed$fname = satmed$fname_new
    # Clean names for plot
    satmed$fname = gsub('wgbs_', '', satmed$fname)
    satmed$fname = gsub('normal_', '', satmed$fname)
    
    # Generate full median dataset fort he line
    nmed = satmed %>% group_by(current_n, fname) %>% summarise(med=median(med), iqr=median(iqr)) %>% data.table()

    ggplot(satmed, aes(x=factor(current_n), y=med, color=fname)) +
    geom_boxplot(position = position_dodge(.5)) +
    #geom_linerange(aes(ymin=med-iqr, ymax=med+iqr),  position=position_dodge(.5)) +
    theme_bw() + geom_point(position=position_dodge(.5)) +
    labs(x="Number of samples in panel", y= "Median percentage\nof DMP CpGs", color="") +
    theme(legend.position='bottom') +
    scale_color_manual(values=c("purple", "orange", 'gray', 'navy', 'red4')) +
    geom_line(aes(x=factor(current_n), y=med, color=fname, group=fname), data=nmed, linetype='dashed', alpha=0.5)
}

plot_val_box <- function(vald){
  vald = vald[vald$DMP!="perc_ns",]
  vald$fname = gsub('wgbs_', '', vald$fname)
  ggplot(vald, aes(x=fname, y=percentage, fill=DMP)) + geom_boxplot(outlier.shape=NA) +
    theme_bw()  + scale_fill_manual(
      labels=c("perc_hyper"="Hyper", "perc_hypo"="Hypo", "perc_ns"="N.S."),
      values=c("Pink", "Skyblue", "white")) +
    #geom_text(aes(label=cg_label, y=max(percentage)+0.5, x=fname), size=3) 
    coord_flip() + labs(x="", y="Percentage CpGs DMP")
  
}



# Write functino to generate boxplot
make_dmp_box <- function(vald){
  ddmp = unique(vald[ , .(n_cg, perc=sum(percentage)), by='item'])
  dboxplot = ggplot(ddmp, aes(x="", y=perc)) +
    theme_bw() +
    geom_boxplot(outlier.shape=NA, fill='red') +
    geom_jitter(alpha=0.5, color="grey", position=position_jitter(width=0.2)) +
    scale_y_continuous(limits=c(0,5)) +
    labs(x="", y="Percentage DMP CpG sites") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
}

# Write function to make dmp fill plot
make_dmp_fill <- function(vald){
dfill = ggplot(vald, aes(x=item, y=percentage, fill=DMP)) +
  geom_col(color='black', position='fill') +
  theme_bw() +
  scale_fill_manual(
    values=c("red", "lightblue")
  ) +
  theme(axis.text.x = element_blank()) +
  labs(x="", y="Fraction of CpGs at DMPs") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
}

# Write function to make density
make_dmp_dens <- function(x, title){ 
dmps = x %>% dplyr::select(contains('mbdiff'))
#dmsamp = dmps
dmsamp = dplyr::sample_n(dmps, min(8e6, nrow(dmps)))
dmlong = tidyr::pivot_longer(dmsamp, cols=1:ncol(dmsamp), names_to='item', values_to='meth')
ddlong = ggplot(dmlong, aes(x=meth, group=item)) + 
  geom_density(aes(fill=item), alpha=0.5, color='purple') +
  theme_bw() +
  labs(x="Methylation Difference\n(Normal - Panel)", y="Density", title=title) +
  theme(legend.position="none") +
  scale_x_continuous(limits=c(-1, 1)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
}

pipe_cplots <- function(x, title="TEST"){
    log_info("Processing {x}")
    # Process data
    dt = qs::qread(x)
    vala = dt[[2]]
    vald = process_val(dt[[1]])
    vald = vald[DMP!="perc_ns",]

    # Combine the three plots together
    cplots <- list(make_dmp_dens(vala, title), make_dmp_box(vald), make_dmp_fill(vald))
    cplots <- patchwork::wrap_plots(cplots, ncol=1)
    # Add title to patchwork and center the title
    cplots <- cplots & theme(plot.title = element_text(hjust = 0.5))#, text=element_text(size=4.5))
    cplots
}

