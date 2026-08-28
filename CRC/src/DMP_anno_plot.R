plot_methylation_info_with_anno <- function(dt, path, sample_id, bulk, return_df=FALSE){

 # Set color code for hyper/hypo 
 cols <- c("gold2", "royalblue2")
 names(cols) <- c("hyper", "hypo")
 
 # Set color code for annotations
 cols2 <- c("black", "grey40", "grey75", "red2", "pink2", "dodgerblue", "forestgreen", 
            "olivedrab3", "lightgreen", "mediumpurple")
 names(cols2) <- c("all_bins", "exon", "intron", "promoter", "prom_CGI", "enhancer","CGI", 
                   "shore", "shelf","repeats")

  tmp <- dt
  if(bulk==FALSE){
       prefix = "purified"
      ylabb = "pure tumour - normal methylation rate"
   } else {
      prefix = "bulk"
      ylabb = "bulk tumour - normal methylation rate"
  }

 # plot fraction hyper/hypo DMRs
 p <- ggplot(data.frame(cols2), aes(x=names(cols2))) + theme_classic() +
      scale_x_discrete(name="", breaks = names(cols2), limits=names(cols2),
                   labels=c("all\nCpG bins", "exon", "intron", "promoter", 
                            "prom\nCGI", "enhancer", "CGI", "shore", "shelf", "repeat")) +
      theme(axis.text.x = element_text(angle = 40, size = 10, vjust=0.75), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
 tmp2 <- data.table(x=names(cols2), y1=0, y2=0, y3=0)
 for(anno in names(cols2)){
   if(anno=="all_bins"){
      flag <- rep(TRUE, nrow(tmp))
   } else if(anno=="prom_CGI"){ 
      flag <- tmp[, promoter=="promoter"&CGI=="CGI"]
   } else {
      nam <- ifelse(anno=="shore", "CGI", 
             ifelse(anno=="shelf", "CGI", 
             ifelse(anno=="repeats", "repeats", anno)))  
       flag <- as.logical(tmp[, .SD, .SDcols=nam]==anno)
   }
   z1 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hyper")/sum(!is.na(DMR_type))]
   z2 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hypo")/sum(!is.na(DMR_type))]
   tmp2[x==anno, c("y1", "y2") := list(z1, z2)]
 
 }

 if(return_df==TRUE){
   return(tmp2)
 } else{
   p1 <- p + scale_y_continuous(name= "hypo-/hypermethylated DMR ratio", breaks = seq(-1,1,0.2), limits=c(-1, 1)) +
           geom_col(data=tmp2, aes(x=x,y=y1,col="hyper", fill="hyper"), alpha=0.5) +
           geom_col(data=tmp2, aes(x=x,y=(-1+y1),col="hypo", fill="hypo"), alpha=0.5) +
           scale_fill_manual(name="anno",values=cols) + 
           scale_color_manual(name="anno",values=cols)
   ggsave(p1, filename = paste0(path, "/", sample_id, "_", prefix, "_rrbs_DMR_summary_plot.pdf"), device="pdf",
             width=6, height = 6, unit="in")
 }
}


plot_methylation_info_with_anno_wgbs <- function(dt, path, sample_id, bulk, return_df=FALSE){

 # Set color code for hyper/hypo 
 cols <- c("gold2", "royalblue2")
 names(cols) <- c("hyper", "hypo")
 
 # Set color code for annotations
 cols2 <- c("black", "grey40", "grey75", "red2", "pink2", "dodgerblue", "forestgreen", 
            "olivedrab3", "lightgreen", "mediumpurple")
 names(cols2) <- c("all_bins", "exon", "intron", "promoter", "prom_CGI", "enhancer","CGI", 
                   "shore", "shelf", "repeat")

  tmp <- dt
  if(bulk==FALSE){
      prefix = "purified"
      ylabb = "pure tumour - normal methylation rate"
   } else {
      prefix = "bulk"
      ylabb = "bulk tumour - normal methylation rate"
  }

 # plot fraction hyper/hypo DMRs
 p <- ggplot(data.frame(cols2), aes(x=names(cols2))) + theme_classic() +
      scale_x_discrete(name="", breaks = names(cols2), limits=names(cols2),
                   labels=c("all\nCpG bins", "exon", "intron", "promoter", 
                            "prom\nCGI", "enhancer", "CGI", "shore", "shelf", "repeat")) +
      theme(axis.text.x = element_text(angle = 40, size = 10, vjust=0.75), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
 tmp2 <- data.table(x=names(cols2), y1=0, y2=0, y3=0)
 for(anno in names(cols2)){
   if(anno=="all_bins"){
      flag <- rep(TRUE, nrow(tmp))
   } else if(anno=="prom_CGI"){ 
      flag <- tmp[, promoter=="promoter" & CGI=="CGI"]
   } else {
      nam <- ifelse(anno=="shore", "CGI", 
             ifelse(anno=="shelf", "CGI",  
             ifelse(anno=="repeat", "repeats", anno)))
       flag <- as.logical(tmp[, .SD, .SDcols=nam]==anno)
   }
   z1 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hyper")/sum(!is.na(DMR_type))]
   z2 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hypo")/sum(!is.na(DMR_type))]
   tmp2[x==anno, c("y1", "y2") := list(z1, z2)]
 }

if(return_df==TRUE){
   return(tmp2)
 } else{
   p1 <- p + scale_y_continuous(name= "hypo-/hypermethylated DMR ratio", breaks = seq(-1,1,0.2), limits=c(-1, 1)) +
           geom_col(data=tmp2, aes(x=x,y=y1,col="hyper", fill="hyper"), alpha=0.5) +
           geom_col(data=tmp2, aes(x=x,y=(-1+y1),col="hypo", fill="hypo"), alpha=0.5) +
           scale_fill_manual(name="anno",values=cols) + 
           scale_color_manual(name="anno",values=cols)
   ggsave(p1, filename = paste0(path, "/", sample_id, "_", prefix, "_wgbs_DMR_summary_plot.pdf"), device="pdf",
             width=6, height = 6, unit="in")
 }
}

plot_DMP_fraction <- function(dt, path, sample_id, bulk, return_df=FALSE){

 # Set color code for hyper/hypo 
 cols <- c("gold2", "royalblue2")
 names(cols) <- c("hyper", "hypo")
 
 # Set color code for annotations
 cols2 <- c("black", "grey40", "grey75", "red2", "pink2", "dodgerblue", "forestgreen", 
            "olivedrab3", "lightgreen", "mediumpurple")
 names(cols2) <- c("all_bins", "exon", "intron", "promoter", "prom_CGI", "enhancer","CGI", 
                   "shore", "shelf","repeats")

  tmp <- dt
  if(bulk==FALSE){
       prefix = "purified"
      ylabb = "pure tumour - normal methylation rate"
   } else {
      prefix = "bulk"
      ylabb = "bulk tumour - normal methylation rate"
  }

 # plot fraction hyper/hypo DMRs
 p <- ggplot(data.frame(cols2), aes(x=names(cols2))) + theme_classic() +
      scale_x_discrete(name="", breaks = names(cols2), limits=names(cols2),
                   labels=c("all\nCpG bins", "exon", "intron", "promoter", 
                            "prom\nCGI", "enhancer", "CGI", "shore", "shelf", "repeat")) +
      theme(axis.text.x = element_text(angle = 40, size = 10, vjust=0.75), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
 tmp2 <- data.table(x=names(cols2), y1=0, y2=0, y3=0)
 for(anno in names(cols2)){
   if(anno=="all_bins"){
      flag <- rep(TRUE, nrow(tmp))
   } else if(anno=="prom_CGI"){ 
      flag <- tmp[, promoter=="promoter"&CGI=="CGI"]
   } else {
      nam <- ifelse(anno=="shore", "CGI", 
             ifelse(anno=="shelf", "CGI", 
             ifelse(anno=="repeats", "repeats", anno)))  
       flag <- as.logical(tmp[, .SD, .SDcols=nam]==anno)
   }
   z1 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hyper")/sum(!is.na(DMR_type) | is.na(DMR_type))]
   z2 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hypo")/sum(!is.na(DMR_type) | is.na(DMR_type))]
   tmp2[x==anno, c("y1", "y2") := list(z1, z2)]
 
 }

 if(return_df==TRUE){
   return(tmp2)
 } else{
   p1 <- p + scale_y_continuous(name= "hypo-/hypermethylated DMR ratio of all CpGs", breaks = seq(-1,1,0.2), limits=c(-1, 1)) +
           geom_col(data=tmp2, aes(x=x,y=y1,col="hyper", fill="hyper"), alpha=0.5) +
           geom_col(data=tmp2, aes(x=x,y=(0-y2),col="hypo", fill="hypo"), alpha=0.5) +
           scale_fill_manual(name="anno",values=cols) + 
           scale_color_manual(name="anno",values=cols)
   ggsave(p1, filename = paste0(path, "/", sample_id, "_", prefix, "_rrbs_DMR_fraction.pdf"), device="pdf",
             width=6, height = 6, unit="in")
 }
}

plot_DMP_fraction_wgbs <- function(dt, path, sample_id, bulk, return_df=FALSE){

 # Set color code for hyper/hypo 
 cols <- c("gold2", "royalblue2")
 names(cols) <- c("hyper", "hypo")
 
 # Set color code for annotations
 cols2 <- c("black", "grey40", "grey75", "red2", "pink2", "dodgerblue", "forestgreen", 
            "olivedrab3", "lightgreen", "mediumpurple")
 names(cols2) <- c("all_bins", "exon", "intron", "promoter", "prom_CGI", "enhancer","CGI", 
                   "shore", "shelf", "repeat")

  tmp <- dt
  if(bulk==FALSE){
      prefix = "purified"
      ylabb = "pure tumour - normal methylation rate"
   } else {
      prefix = "bulk"
      ylabb = "bulk tumour - normal methylation rate"
  }

 # plot fraction hyper/hypo DMRs
 p <- ggplot(data.frame(cols2), aes(x=names(cols2))) + theme_classic() +
      scale_x_discrete(name="", breaks = names(cols2), limits=names(cols2),
                   labels=c("all\nCpG bins", "exon", "intron", "promoter", 
                            "prom\nCGI", "enhancer", "CGI", "shore", "shelf", "repeat")) +
      theme(axis.text.x = element_text(angle = 40, size = 10, vjust=0.75), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
 tmp2 <- data.table(x=names(cols2), y1=0, y2=0, y3=0)
 for(anno in names(cols2)){
   if(anno=="all_bins"){
      flag <- rep(TRUE, nrow(tmp))
   } else if(anno=="prom_CGI"){ 
      flag <- tmp[, promoter=="promoter"&CGI=="CGI"]
   } else {
      nam <- ifelse(anno=="shore", "CGI", 
             ifelse(anno=="shelf", "CGI",  
             ifelse(anno=="repeat", "repeats", anno)))  
       flag <- as.logical(tmp[, .SD, .SDcols=nam]==anno)
   }
   z1 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hyper")/sum(!is.na(DMR_type) | is.na(DMR_type))]
   z2 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hypo")/sum(!is.na(DMR_type) | is.na(DMR_type))]
   tmp2[x==anno, c("y1", "y2") := list(z1, z2)]
 
 }

if(return_df==TRUE){
   return(tmp2)
   } else {
   p1 <- p + scale_y_continuous(name= "hypo-/hypermethylated DMR ratio of all CpGs", breaks = seq(-1,1,0.2), limits=c(-1, 1)) +
           geom_col(data=tmp2, aes(x=x,y=y1,col="hyper", fill="hyper"), alpha=0.5) +
           geom_col(data=tmp2, aes(x=x,y=(0-y2),col="hypo", fill="hypo"), alpha=0.5) +
           scale_fill_manual(name="anno",values=cols) + 
           scale_color_manual(name="anno",values=cols)
   ggsave(p1, filename = paste0(path, "/", sample_id, "_", prefix, "_wgbs_DMR_fraction.pdf"), device="pdf",
             width=6, height = 6, unit="in")
   }
}