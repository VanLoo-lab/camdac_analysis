#load libs
x <- c("DescTools", "eulerr", "VennDiagram", "gridExtra", "ggplot2", 
       "data.table", "fst", "limma", "dplyr","ggforce")
# x <- c(x, "umap", "GenomicRanges")
invisible(lapply(x, function(y) 
  suppressWarnings(suppressMessages(library(y,character.only = TRUE,
                                            quietly=TRUE,warn.conflicts = FALSE)))))
rm(x)

# # Set vars 
# patient_ids = c("LTX092", "LTX064", "LTX058") 
# effect_size=0.2
# prob=0.01
# path_to_CAMDAC="~/"
# control_id="N"
# min_DMP_counts_in_DMR=5
# min_consec_DMP_in_DMR=3
# n_cores = 12

parent_dir <- "~/rdeocd/txMeth/"
out_dir <- "~/data/CAMDAC_paper/figure3DE/outputs"

# load patient metadata
f_nm <- paste0(parent_dir, 
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt")
patient_info <- fread(f_nm, sep = "\t")

# set run params
reseg <- TRUE
effect_size <- 0.2

# for (patient_id1 in patient_ids) {
for (patient_id1 in "LTX064") {
  path1 <- paste0("~/rdeocd/txMeth/",patient_id1, "/")
  sample_ids1 <- unique(patient_info[Patient == patient_id1, Sample])
  
  cols <- c("chromosome", "start", "end", "m_t_low", "m_t_high", "m_t_corr",
            "m_b_low", "m_b_high", "m_b"#, "m_n_low", "m_n_high", "m_n",
            )#"nA", "nB", "CN"

  df_list1 <- NULL
  for (i in seq_along(sample_ids1)) {
    sample_id <- sample_ids1[i]
    f_nm <- paste0(path1, "/Methylation/",sample_id,"/re_segmented/purified_tumour_DMP_calls.RData")
    if (reseg == FALSE) {
      f_nm <- paste0(path1, "/Methylation/",sample_id,"/purified_tumour_DMP_calls.RData")
    }
    load(f_nm)
    df_mt <- data.table(df_mt[,cols])
    colnames(df_mt)[4:length(cols)] <- paste(cols[4:length(cols)],patient_id1,sample_id,sep="_")
    df_list1[[length(df_list1)+1]] <- df_mt
    if(length(df_list1)>1){
      df1 <- Reduce(function(...) merge(..., by=c("chromosome","start","end"),all = T), df_list1) 
      df_list1 <- list(df1)
      rm(df1)
    }
    rm(df_mt)
  }
  df1 <- df_list1
  rm(df_list1)
  
  for (patient_id2 in "LTX064") { #c("LTX064", "LTX038")
    path2 <- paste0("~/rdeocd/txMeth/", patient_id2, "/")
    sample_ids2 <- unique(patient_info[Patient == patient_id2, Sample])
    
    if(patient_id1 == patient_id2){
      df2 <- df1
    } else {
    df_list2 <- NULL
    
    # Load all tumour samples for patient 2 with DMP annotated
    for (i in seq_along(sample_ids2)) {
      sample_id <- sample_ids2[i]
      f_nm <- paste0(path2, "/Methylation/",sample_id,"/re_segmented/purified_tumour_DMP_calls.RData")
      if(reseg==FALSE){f_nm <- paste0(path2, "/Methylation/",sample_id,"/purified_tumour_DMP_calls.RData")}
      load(f_nm)
      df_mt <- data.table(df_mt[,cols])
      colnames(df_mt)[4:length(cols)] <- paste(cols[4:length(cols)],patient_id2,sample_id,sep="_")
      df_list2[[length(df_list2)+1]] <- df_mt
      if(length(df_list2)>1){
        df2 <- Reduce(function(...) merge(..., by=c("chromosome","start","end"),all = T), df_list2) 
        df_list2 <- list(df2)
        rm(df2)
      }
      rm(df_mt)
    }
    df2 <- df_list2
    rm(cols, f_nm, sample_id, i)
    }

    df <- Reduce(function(...) merge(..., by=c("chromosome","start","end"),all = T), list(df1, df2))
    dt <- data.table(df)
    DMPs_mt <- as.data.table(df[, 1:3])
    rm(df1, df2)

    for(sample_id1 in sample_ids1){
      for(sample_id2 in sample_ids2){
        if( patient_id1 != patient_id2 |
          (patient_id1 == patient_id2 & which(sample_ids2 == sample_id2) >
            which(sample_ids1==sample_id1))) {
        cat(patient_id1, "-", sample_id1, " vs. ",
            patient_id2, "-", sample_id2, "\n")
        raw <- rep(NA, nrow(DMPs_mt))
        cols1 <- paste(c("m_t_low", "m_t_high"),patient_id1, sample_id1, sep="_")
        cols2 <- paste(c("m_t_low", "m_t_high"),patient_id2, sample_id2, sep="_")

        if (patient_id1 == patient_id2) {
          cols1 <- paste(cols1, "x", sep = ".")
          cols2 <- paste(cols2, "y", sep = ".")
        }
        tmp <- rowSums(dt[, is.na(.SD), .SDcols = c(cols1[1], cols2[1])])
        idx <- c(1:nrow(df))[tmp == 0]
        rm(tmp)
        
        raw[idx] <- Overlap(df[idx, cols1], df[idx, cols2]) == 0
        if (patient_id1 != patient_id2) {
          diff <- df[,paste(c("m_t_corr"), patient_id1, sample_id1, sep="_")] -
            df[,paste(c("m_t_corr"), patient_id2, sample_id2, sep = "_")]
        } else {
          diff <- df[,paste0("m_t_corr_", patient_id1, "_", sample_id1, ".x")] -
            df[,paste0("m_t_corr_", patient_id2, "_", sample_id2, ".y")]
        }
        corr <- abs(diff) >= effect_size & raw == TRUE

        # cols <- c(paste("raw",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"),
        #           paste("corrected",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"))
        cols <- c(paste("corrected",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"))
        DMPs_mt[, as.character(cols) := list(NA)]
        # DMPs_mt[,cols] <- cbind(raw, corr)
        DMPs_mt[, as.character(cols) := as.list(corr)]
        }
      }
    }
    rm(raw, diff, corr)

    DMPs_mb <- as.data.table(df[,1:3])
    for(sample_id1 in sample_ids1){
      for(sample_id2 in sample_ids2){
        if(patient_id1!=patient_id2|(patient_id1==patient_id2&which(sample_ids2==sample_id2)>which(sample_ids1==sample_id1))){
        cat(patient_id1, "-", sample_id1, " vs. ",
            patient_id2, "-", sample_id2, "\n")
        raw <- rep(NA, nrow(DMPs_mb))  
        cols1 <- paste(c("m_b_low", "m_b_high"),patient_id1, sample_id1, sep="_")
        cols2 <- paste(c("m_b_low", "m_b_high"),patient_id2, sample_id2, sep="_")
        if(patient_id1==patient_id2){
          cols1 <- paste(cols1, "x", sep=".")
          cols2 <- paste(cols2, "y", sep=".")
        }
        
        tmp <- rowSums(dt[,is.na(.SD), .SDcols=c(cols1[1],cols2[1])])
        idx <- c(1:nrow(df))[tmp==0] 
        rm(tmp)
        
        raw[idx] <- Overlap(df[idx,cols1],df[idx,cols2])==0
        if(patient_id1!=patient_id2){
          diff <- df[,paste(c("m_b"),patient_id1, sample_id1, sep="_")] - df[,paste(c("m_b"),patient_id2, sample_id2, sep="_")]
        } else {
          diff <- df[,paste0("m_b_",patient_id1,"_", sample_id1, ".x")] - df[,paste0("m_b_",patient_id2,"_", sample_id2, ".y")]
        }
        corr <- abs(diff) >= effect_size & raw == T
        
        # cols <- c(paste("raw",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"), 
        #           paste("corrected",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"))
        cols <- c(paste("corrected",patient_id1,sample_id1,patient_id2,sample_id2, sep="_"))
        DMPs_mb[, as.character(cols) := NA]
        # DMPs_mb[,cols] <- cbind(raw, corr)
        DMPs_mb[, as.character(cols) := as.list(corr)]
        }
      }
    }
    rm(raw, diff, corr)

    # get all DMPs loci DMPs between at least one pair in bulk or pure tumour-tumour DMP set
    DMPs <- na.omit(cbind.data.frame(DMPs_mt, DMPs_mb))
    cols_corrected <- grep(("corrected_"),colnames(DMPs_mt))
    # cols_raw <- grep(("raw_"),colnames(DMPs_mt))
    # 
    # df_plot_corr <- list(); df_plot_raw <- list()
    # df_plot_corr_rates <- list(); df_plot_raw_rates <- list()
    venn_corr <- list()
    #venn_raw <- list()
    # 
    max_DMPs_corr <- max(DMPs_mt[,lapply(.SD, function(x) sum(as.numeric(x), na.rm=T)), .SDcols=cols_corrected],
                         DMPs_mb[,lapply(.SD, function(x) sum(as.numeric(x), na.rm=T)), .SDcols=cols_corrected])
    # max_DMPs_raw <- max(colSums(DMPs_mt[,cols_raw], na.rm=T),colSums(DMPs_mb[,cols_raw], na.rm=T)) 

    venn_overlap_list <- list()
    for (sample_id1 in sample_ids1) {
      for (sample_id2 in sample_ids2) {
        if (patient_id1 != patient_id2 |
           (patient_id1 == patient_id2 &
            which(sample_ids2 == sample_id2) > which(sample_ids1 == sample_id1))) {
        cols <- c(grep(eval(paste("corrected", patient_id1, sample_id1, patient_id2, sample_id2, 
                                  sep="_")),colnames(DMPs_mt)))
        tmp <- na.omit(data.frame(mt = DMPs_mt[, unlist(.SD), .SDcols = c(cols)],
                                  mb = DMPs_mb[, unlist(.SD), .SDcols = c(cols)]))
        rownames(tmp) <- 1:nrow(tmp)
        
        # get counts ffor venn diagram
        tmp2 <- euler(tmp, shape = "ellipse")
        venn_overlap_list[[length(venn_overlap_list) + 1]] <-
          data.table(sample_id1 = sample_id1, sample_id2 = sample_id2,
                     perc_mt_also_mb = tmp2$original[3] / sum(tmp2$original[c(1,3)]),
                     perc_mb_also_mt = tmp2$original[3] / sum(tmp2$original[c(2,3)]),
                     mt_only = tmp2$original[1],
                     mb_only = tmp2$original[2],
                     both = tmp2$original[3])
        
        # plot venn counts
        f_prefix <- paste0(patient_id1, "-", sample_id1, "_", patient_id2, "-", sample_id2)
        f_nm <- paste0(f_prefix, "_mt_mt_vs_mb_mb_DMPs.pdf")
        pdf(file.path(out_dir, f_nm), width = 3, height = 2.5)
        print(plot(tmp2, quantities = TRUE))
        dev.off()
        Sys.sleep(1)
        rm(tmp2)
        
        # # estract CN difference
        # cols_tmp <- paste0("CN_", c(sample_id1,sample_id2))
        # tmp <- inner_join(tmp,df[,c(colnames(DMPs_mt[,1:3]),cols_tmp)])
        # tmp$CN_diff <- ifelse(tmp[,8]/tmp[,9] > 1, tmp[,8]/tmp[,9], tmp[,9]/tmp[,8])
        # tmp$CN_diff <- cut(tmp$CN_diff, c(0,1,2,Inf), labels=c("A","B","C"))
        # tmp$CN_diff <- abs(tmp[,8] - tmp[,9])
        # tmp$CN_diff <- ifelse(tmp$CN_diff >=2, "≥2", tmp$CN_diff)
        # tmp[,1:3] <- NULL
        # rm(cols_tmp)

        # # Get venn counts for all 3 CN diff categories (0,1,2)
        # labs <- factor(sort(unique(tmp$CN_diff)), levels=c("A","B","C"))
        # CN_results[[counter]] <- data.frame(label= labs,Counts_FP=numeric(length = 3), Counts_TP=numeric(length = 3),Counts_FN=numeric(length = 3))
        # for (n in labs){
        #   tmp_tmp <- tmp[tmp$CN_diff==n,]
        #   vdc <- vennCounts(tmp_tmp[,c(2,4)])
        #   class(vdc) <- 'matrix'
        #   df.vdc <- as.data.frame(vdc)[-1,]
        #   df.vdc.ordered <- df.vdc[c(1,3,2),]
        # 
        #   # label venn counts
        #   df.vdc.ordered$label = rowSums(df.vdc.ordered[,1:2])
        #   df.vdc.ordered$label = ifelse(df.vdc.ordered[,grepl("_mt", colnames(df.vdc.ordered))]==1&df.vdc.ordered$label==1, "false negatives",
        #                                 ifelse(df.vdc.ordered$label==2, "true positives","false positives"))
        #   # Store results
        #   CN_results[[counter]][which(sort(unique(tmp$CN_diff))==n),1] <- n
        #   CN_results[[counter]][which(sort(unique(tmp$CN_diff))==n),2:4] <- c(df.vdc.ordered[df.vdc.ordered$label=="false positives","Counts"],
        #                                                                       df.vdc.ordered[df.vdc.ordered$label=="true positives","Counts"],
        #                                                                       df.vdc.ordered[df.vdc.ordered$label=="false negatives","Counts"])
        #   };rm(tmp_tmp)

        # Get venn counts (corrected calls)
        df.venn <- data.frame(x = c(0, NA),y = c(1, 1),labels = c( 'mb1/mb2\n DMPs','mt1/mt2\n DMPs'))
        vdc <- vennCounts(tmp[,c(1, 2)])
        class(vdc) <- 'matrix'
        df.vdc <- as.data.frame(vdc)[-1,]
        df.vdc.ordered <- df.vdc[c(1,3,2),] %>%
            mutate(x = NA, y = c(1, 1, 1))
        rm(vdc,df.vdc)

        # label venn counts
        df.vdc.ordered$label = rowSums(df.vdc.ordered[, 1:2])
        df.vdc.ordered$label = ifelse(df.vdc.ordered[,grepl("mt", colnames(df.vdc.ordered))]==1&df.vdc.ordered$label==1, "false negatives",
                                      ifelse(df.vdc.ordered$label==2, "true positives","false positives"))

        # get ven circle area setting the mt1-mt2 circle area (circle B) to pi (i.e. radius = 1)
        coordinates <- data.frame(label= c("mt","shared","mb"),
                                  Counts = c(sum(df.vdc.ordered[df.vdc.ordered$label!="false positives","Counts"]),
                                            df.vdc.ordered[df.vdc.ordered$label=="true positives","Counts"],
                                            sum(df.vdc.ordered[df.vdc.ordered$label!="false negatives","Counts"])), stringsAsFactors = F)
        coordinates$area <- c(pi*(coordinates$Counts[1]/max_DMPs_corr), pi*(coordinates$Counts[2]/max_DMPs_corr), pi*(coordinates$Counts[3]/max_DMPs_corr))
        coordinates$rad <- sqrt(coordinates$area/pi); coordinates$rad[2] <- NA
        coordinates$rad2 <- coordinates$area/pi ; coordinates$rad2[2] <- NA
        df.venn$rad <- c(coordinates$rad[3],coordinates$rad[1])

        if(coordinates$area[2]>0.5*coordinates$area[1]){
          s <- seq(-min(coordinates$rad, na.rm = T)+1E-12,min(coordinates$rad, na.rm = T)-1E-12,0.0001)
          tmp2 <- pi*coordinates$rad2[1] - (coordinates$rad2[1]*asin(s/coordinates$rad[1])) + s*sqrt(coordinates$rad2[1]-(s^2)) +
              (coordinates$rad2[3]*asin(s/coordinates$rad[3])) - s*sqrt(coordinates$rad2[3]-(s^2)) - coordinates$area[2]
          y_coord <- s[which(tmp2==min(tmp2[tmp2>=0]))]
          rm(tmp2, s)
          # fct <- function(x, y, coordinates) {
          #   pi*coordinates$rad2[1] - (coordinates$rad2[1]*asin(x/coordinates$rad[1])) + x*sqrt(coordinates$rad2[1]-(x^2)) +
          #     (coordinates$rad2[3]*asin(x/coordinates$rad[3])) - x*sqrt(coordinates$rad2[3]-(x^2)) - y
          # }
        }
        if(coordinates$area[2]<=0.5*coordinates$area[1]){
          s <- seq(-min(coordinates$rad, na.rm = T)+1E-12,min(coordinates$rad, na.rm = T)-1E-12,0.0001)
          tmp2 <- (coordinates$rad2[1]*asin(s/coordinates$rad[1])) - s*sqrt(coordinates$rad2[1]-(s^2)) +
              (coordinates$rad2[3]*asin(s/coordinates$rad[3])) - s*sqrt(coordinates$rad2[3]-(s^2)) - coordinates$area[2]
          y_coord <- s[which(tmp2==min(tmp2[tmp2>=0]))]
          rm(tmp2, s)
          # fct <- function(x, y, coordinates) {
          #   (coordinates$rad2[1]*asin(x/coordinates$rad[1])) + (coordinates$rad2[3]*asin(x/coordinates$rad[3])) -
          #     x*sqrt(coordinates$rad2[1]-(x^2)) - x*sqrt(coordinates$rad2[3]-(x^2)) - y
          # }
        }

        fct2 <- function(x, x_coord, radx, rady) { x^2 - 2*x_coord*x + radx - rady }

        # find the x and y coordinate of the intersecting points giving circle A center aligned with (0,0) and circle B with (0,d)
        # where d is the distance between circle A and circle B
        # y_coord <- abs(uniroot(fct,y=coordinates$area[2],lower=-min(coordinates$rad, na.rm = T)+0.0001,upper=min(coordinates$rad, na.rm = T)-0.0001, f.upper = coordinates$area[2])$root)
        x_coord <- sqrt(coordinates$rad2[3] - (y_coord^2))

        # find the distance between the two circles
        radx <- coordinates$rad2[3]
        rady <-coordinates$rad2[1]
        d <- UnirootAll(fct2,x_coord=x_coord, radx=radx, rady=rady, lower=-min(coordinates$rad, na.rm = T)+1E-6,upper=sum(coordinates$rad, na.rm = T))
        if(paste0(sample_id1, "_", sample_id2)%in%c(
          "R1_R5", "R1_R7", "R2_R3", "R2_R4", "R2_R7", "R3_R5", "R4_R5", "R3_R7", "R4_R7", "R5_R7")){ d <- d[1]}
        # for R1 and R3 there are 2 roots.. 0.6157418 and 0.6284104, use 0.622
        # R1/R4 0.1593 R5/R5 0.2125 R5/R6 0.3211 R6/R1 0.5950
        if(length(d)==0){
          idx <-  which.min(abs(fct2(seq(0,1,0.0001), x_coord, radx, rady)))
          d <- seq(0,1,0.0001)[idx] ; rm(idx)
        } else {
          #d <- ifelse(coordinates$Counts[coordinates$label=="shared"]>df.vdc.ordered$Counts[3], min(d), max(d)) 
          d <- d[which.min(abs(coordinates$rad[1]+sqrt(coordinates$area/pi)[2]-d))]
        }
        df.venn$x[2] <- d
        stopifnot(is.finite(d))

        df.vdc.ordered$x <- c(0-(0.5*coordinates$rad[3]),x_coord,d+((2*coordinates$rad[1]-d)/2))
        df.vdc.ordered$y <- c(1+y_coord,1,1+y_coord)

        # Plot venn diagram
        out_file <- file.path(out_dir,
          paste(sample_id1, sample_id2, "mt_mt_vs_mb_mb_DMPs_corrected.pdf", sep = "_")
          )
        ggplot(df.venn) + coord_fixed() + theme_void() + theme(legend.position = 'right') +
          geom_circle(aes(x0 = x, y0 = y, r = rad, color = labels, fill=labels), alpha = .1, size=2) +
          scale_colour_manual(name="legend", values = c("lightgrey", 'dodgerblue')) + 
          scale_fill_manual(name="legend",values = c("lightgrey", 'dodgerblue'))+
          labs(fill = NULL) + 
          annotate("text", x = df.vdc.ordered$x, y = df.vdc.ordered$y, label = df.vdc.ordered$Counts, size = 5)
        ggsave(file = out_file, device = "pdf", width=4.5, height = 2.5, units = "in")

        # save venn info
        venn_corr[[length(venn_corr)+1]] <- df.venn
        venn_corr[[length(venn_corr)]]$sample_id1 <- sample_id1
        venn_corr[[length(venn_corr)]]$sample_id2 <- sample_id2

        # # save venn counts
        # df_plot_corr[[counter]] <- df.vdc.ordered[,c("Counts", "label", "x", "y")]
        # df_plot_corr[[counter]]$sample_id1 <- sample_id1 ; df_plot_corr[[counter]]$sample_id2 <- sample_id2
        # df_plot_corr[[counter]]$sum <- sum(df_plot_corr[[counter]]$Counts)
        # df_plot_corr[[counter]]$prop <- df_plot_corr[[counter]]$Counts/df_plot_corr[[counter]]$sum
        # df_plot_corr[[counter]]$prop_min <- c(0, cumsum(df_plot_corr[[counter]]$prop[1:2]))
        # df_plot_corr[[counter]]$prop_max <- c(cumsum(df_plot_corr[[counter]]$prop[1:2]),1)
        # 
        # # TP, TN, FN, FP
        # df_plot_corr_rates[[counter]] <- df_plot_corr[[counter]][,c("Counts", "label")]
        # df_plot_corr_rates[[counter]] <- rbind(df_plot_corr_rates[[counter]], cbind.data.frame(Counts=nrow(DMPs)-sum(df_plot_corr_rates[[counter]]$Counts), label="true negatives"))
        # 
        # # Get venn counts (corrected calls)
        # df.venn <- data.frame(x = c(0, NA),y = c(1, 1),labels = c( 'mb1/mb2\n DMPs','mt1/mt2\n DMPs'))
        # vdc <- vennCounts(tmp[,c(1,3)])
        # class(vdc) <- 'matrix'
        # df.vdc <- as.data.frame(vdc)[-1,] 
        # df.vdc.ordered <- df.vdc[c(1,3,2),]%>%
        #   mutate(x = NA, y = c(1, 1, 1))
        # rm(vdc,df.vdc)
        # 
        # # label venn counts
        # df.vdc.ordered$label = rowSums(df.vdc.ordered[,1:2])
        # df.vdc.ordered$label = ifelse(df.vdc.ordered[,grepl("_mt", colnames(df.vdc.ordered))]==1&df.vdc.ordered$label==1, "false negatives",
        #                               ifelse(df.vdc.ordered$label==2, "true positives","false positives"))
        # 
        # # get ven circle area setting the mt1-mt2 circle area (circle B) to pi (i.e. radius = 1) 
        # coordinates <- data.frame(label= c("mt","shared","mb"), 
        #                           Counts = c(sum(df.vdc.ordered[df.vdc.ordered$label!="false positives","Counts"]),
        #                                      df.vdc.ordered[df.vdc.ordered$label=="true positives","Counts"],
        #                                      sum(df.vdc.ordered[df.vdc.ordered$label!="false negatives","Counts"])), stringsAsFactors = F)
        # coordinates$area <- c(pi*(coordinates$Counts[1]/max_DMPs_corr), pi*(coordinates$Counts[2]/max_DMPs_corr), pi*(coordinates$Counts[3]/max_DMPs_corr))   
        # coordinates$rad <- sqrt(coordinates$area/pi); coordinates$rad[2] <- NA
        # coordinates$rad2 <- coordinates$area/pi ; coordinates$rad2[2] <- NA
        # df.venn$rad <- c(coordinates$rad[3],coordinates$rad[1])
        # 
        # if(coordinates$area[2]>0.5*coordinates$area[1]){
        #   s <- seq(-min(coordinates$rad, na.rm = T)+1E-12,min(coordinates$rad, na.rm = T)-1E-12,0.0001)
        #   tmp2 <- pi*coordinates$rad2[1] - (coordinates$rad2[1]*asin(s/coordinates$rad[1])) + s*sqrt(coordinates$rad2[1]-(s^2)) +
        #     (coordinates$rad2[3]*asin(s/coordinates$rad[3])) - s*sqrt(coordinates$rad2[3]-(s^2)) - coordinates$area[2]
        #   y_coord <- s[which(tmp2==min(tmp2[tmp2>=0]))]
        #   rm(tmp2, s)
        #   # fct <- function(x, y, coordinates) { 
        #   #   pi*coordinates$rad2[1] - (coordinates$rad2[1]*asin(x/coordinates$rad[1])) + x*sqrt(coordinates$rad2[1]-(x^2)) + 
        #   #     (coordinates$rad2[3]*asin(x/coordinates$rad[3])) - x*sqrt(coordinates$rad2[3]-(x^2)) - y 
        #   # }
        # }
        # if(coordinates$area[2]<=0.5*coordinates$area[1]){
        #   s <- seq(-min(coordinates$rad, na.rm = T)+1E-12,min(coordinates$rad, na.rm = T)-1E-12,0.0001)
        #   tmp2 <- (coordinates$rad2[1]*asin(s/coordinates$rad[1])) - s*sqrt(coordinates$rad2[1]-(s^2)) +
        #     (coordinates$rad2[3]*asin(s/coordinates$rad[3])) - s*sqrt(coordinates$rad2[3]-(s^2)) - coordinates$area[2]
        #   y_coord <- s[which(tmp2==min(tmp2[tmp2>=0]))]
        #   rm(tmp2, s)
        #   # fct <- function(x, y, coordinates) { 
        #   #   (coordinates$rad2[1]*asin(x/coordinates$rad[1])) + (coordinates$rad2[3]*asin(x/coordinates$rad[3])) -
        #   #     x*sqrt(coordinates$rad2[1]-(x^2)) - x*sqrt(coordinates$rad2[3]-(x^2)) - y 
        #   # }
        # }   
        # 
        # 
        # fct2 <- function(x, x_coord, radx, rady) { x^2 - 2*x_coord*x + radx - rady }
        # # seq(-sum(coordinates$rad, na.rm = T),sum(coordinates$rad, na.rm = T),0.0001)^2 + 2*x_coord*seq(-sum(coordinates$rad, na.rm = T),
        # #      sum(coordinates$rad, na.rm = T),0.01) + coordinates$rad2[3] - coordinates$rad2[1]
        # 
        # # find the x and y coordinate of the intersecting points giving circle A center aligned with (0,0) and circle B with (0,d) 
        # # where d is the distance between circle A and circle B
        # #y_coord <- abs(uniroot(fct,y=coordinates$area[2],coordinates=coordinates,lower=-min(coordinates$rad, na.rm = T)+1E-12,upper=min(coordinates$rad, na.rm = T)-1E-12, f.upper = coordinates$area[2])$root)
        # x_coord <- sqrt(coordinates$rad2[3] - (y_coord^2))
        # 
        # # find the distance between the two circles
        # radx <- coordinates$rad2[3]
        # rady <-coordinates$rad2[1]
        # d <- UnirootAll(fct2,x_coord=x_coord, radx=radx, rady=rady, lower=-min(coordinates$rad, na.rm = T)+1E-6,upper=sum(coordinates$rad, na.rm = T))
        # d <- ifelse(coordinates$Counts[coordinates$label=="shared"]>df.vdc.ordered$Counts[3], min(d), max(d)) 
        # stopifnot(is.finite(d))
        # df.venn$x[2] <- d
        # 
        # if(coordinates$area[2]<=0.5*coordinates$area[1]){df.vdc.ordered$x <-c(0-(0.5*coordinates$rad[3]),x_coord,d+((2*coordinates$rad[1]-d)/2))}
        # if(coordinates$area[2]>0.5*coordinates$area[1]){df.vdc.ordered$x <- c(0-(0.5*coordinates$rad[3]),d,x_coord+0.1*coordinates$rad[1])}
        # df.vdc.ordered$y <- c(1+y_coord,1,1+y_coord)
        # 
        # # Plot venn diagram
        # ggplot(df.venn) + coord_fixed() + theme_void() + theme(legend.position = 'right') +
        #   geom_circle(aes(x0 = x, y0 = y, r = rad, color = labels, fill=labels), alpha = .1, size=2) +
        #   scale_colour_manual(name="legend",values = c("grey20", 'darkblue')) + scale_fill_manual(name="legend",values = c("grey20", 'darkblue'))+
        #   labs(fill = NULL) + annotate("text", x = df.vdc.ordered$x, y = df.vdc.ordered$y, label = df.vdc.ordered$Counts, size = 5)+
        #   ggsave(file = paste0(out_dir,sample_id1,"_",sample_id2,"_mt_mt_vs_mb_mb_DMPs_raw.pdf"),device = "pdf", width=4.5, height = 2.5, units = "in")
        # 
        # # save venn info
        # venn_raw[[counter]] <- df.venn
        # venn_raw[[counter]]$sample_id1 <- sample_id1 ; venn_raw[[counter]]$sample_id2 <- sample_id2
        # 
        # # save venn counts
        # df_plot_raw[[counter]] <- df.vdc.ordered[,c("Counts", "label", "x", "y")]
        # df_plot_raw[[counter]]$sample_id1 <- sample_id1 ; df_plot_raw[[counter]]$sample_id2 <- sample_id2
        # df_plot_raw[[counter]]$sum <- sum(df_plot_raw[[counter]]$Counts)
        # df_plot_raw[[counter]]$prop <- df_plot_raw[[counter]]$Counts/df_plot_raw[[counter]]$sum
        # df_plot_raw[[counter]]$prop_min <- c(0, cumsum(df_plot_raw[[counter]]$prop[1:2]))
        # df_plot_raw[[counter]]$prop_max <- c(cumsum(df_plot_raw[[counter]]$prop[1:2]),1)
        # 
        # # TP, TN, FN, FP
        # df_plot_raw_rates[[counter]] <- df_plot_raw[[counter]][,c("Counts", "label")]
        # df_plot_raw_rates[[counter]] <- rbind(df_plot_raw_rates[[counter]], cbind.data.frame(Counts=nrow(DMPs)-sum(df_plot_raw_rates[[counter]]$Counts), label="true negatives"))
        }
      }
    }
  }
}
rm(df.venn, df.vdc.ordered, counter, tmp, cols, coordinates, x_coord, y_coord, d)

df_overlaps <- rbindlist(venn_overlap_list) 
median(df_overlaps$perc_mt_also_mb)
median(df_overlaps$perc_mb_also_mt)
save(df_overlaps, file = file.path(out_dir, "intra_tumour_DMPs.RData"))
rm(venn_overlap_list)

# merge dfs for plotting
df_plot <-  bind_rows(venn_corr)
samples <- unique(paste(df_plot$sample_id1,df_plot$sample_id2,sep="_"))

rho1 <- patient_info[Patient==patient_id1, .(Sample, rho_reseg)]
rho2 <- patient_info[Patient==patient_id2, .(Sample, rho_reseg)]
df_plot$sample_id1 <- factor(df_plot$sample_id1, levels=rho1$Sample[order(rho1$rho_reseg)], ordered = TRUE)
df_plot$sample_id2 <- factor(df_plot$sample_id2, levels=rho2$Sample[order(rho2$rho_reseg)], ordered = TRUE)

# Plot summary venn diagram
ggplot(data=df_plot[is.finite(df_plot$x)&!grepl("raw", df_plot$labels),]) + coord_fixed() + theme_void() +
  facet_wrap(sample_id1~sample_id2, ncol = length(sample_ids2), nrow=length(sample_ids1), drop = F) +
  geom_circle(aes(x0 = x, y0 = y, r = rad, color = labels, fill = labels), alpha = .25, size=.5) +
  theme(legend.position = 'right',axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(),legend.text=element_text(size=8),
        axis.ticks = element_blank(), strip.background = element_blank(),
        #strip.text.x = element_blank()
  ) +
  scale_colour_manual(name="legend", values = c("mb1/mb2\n DMPs"="lightgrey","mt1/mt2\n DMPs"="dodgerblue",
                                                "mb1/mb2\n raw DMPs"="grey20","mt1/mt2\n raw DMPs"="darkblue")) +
  scale_fill_manual(name="legend",values = c("mb1/mb2\n DMPs"="lightgrey","mb1/mb2\n raw DMPs"="grey20",
                                             "mt1/mt2\n DMPs"="dodgerblue","mt1/mt2\n raw DMPs"="darkblue")) +
  ggsave(file = paste0(out_dir, "/summary_figure_venn_corrected.pdf"),device = "pdf", width=5, height = 6.5, units = "in")


# # merge dfs for plotting
# df_plot <-  bind_rows(df_plot_corr) ; tmp <- bind_rows(df_plot_raw)
# df_plot <- rbind(df_plot,cbind(tmp[,1:4], sample_id1=tmp$sample_id2, sample_id2=tmp$sample_id1, tmp[,7:10]))
# 
# # set pie area relative to number of DMPs
# #df_plot$pie_area <- sqrt(df_plot$sum/max(df_plot$sum))/pi
# df_plot$pie_area <- pi*(df_plot$sum/max(df_plot$sum))
# 
# # order facet
# df_plot$sample_id1 <- factor(df_plot$sample_id1, levels = sample_ids, ordered=T)
# df_plot$sample_id2 <- factor(df_plot$sample_id2, levels = sample_ids, ordered=T)
# 
# # modify labels for raw DMPs
# df_plot$label[((nrow(df_plot)/2)+1):nrow(df_plot)] <- paste0(df_plot$label[((nrow(df_plot)/2)+1):nrow(df_plot)], " raw")
# 
# # plot summary pies
# ggplot(df_plot, aes(fill=label, color=label))+ theme_classic() + facet_wrap(sample_id1~sample_id2,nrow = length(sample_ids), ncol=length(sample_ids), drop = F) +
#   geom_rect(aes(xmin=prop_min,xmax=prop_max,ymin=0,ymax=pie_area),alpha=0.25,stat = "identity",na.rm = T) +   coord_polar("x") +
#   theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),legend.text=element_text(size=8), axis.ticks = element_blank(), strip.background = element_blank(),
#         strip.text.x = element_blank()) +
#   scale_fill_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                           "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   scale_color_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                            "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   ggsave(file = paste0(path, "Methylation/mb_mb_vs_mt_mt_DMPs/summary_figure.pdf"),device = "pdf", width=5, height = 5.5, units = "in")
# 
# # get FN/TP+FN and FP/TP+FP
# prop1 <- numeric(length(df_plot_corr_rates))
# prop2 <- numeric(length(df_plot_corr_rates))
# for(i in 1:length(df_plot_corr_rates)){
#   prop1[i] <- df_plot_corr_rates[[i]]$Counts[1]/(df_plot_corr_rates[[i]]$Counts[1]+df_plot_corr_rates[[i]]$Counts[2])
#   prop2[i] <- df_plot_corr_rates[[i]]$Counts[3]/(df_plot_corr_rates[[i]]$Counts[1]+df_plot_corr_rates[[i]]$Counts[3])
# }
# 
# # get FN/TP+FN and FP/TP+FP
# for(i in 1:length(df_plot_raw_rates)){
#   prop1[length(df_plot_raw_rates)+i] <- df_plot_raw_rates[[i]]$Counts[1]/(df_plot_raw_rates[[i]]$Counts[1]+df_plot_raw_rates[[i]]$Counts[2])
#   prop2[length(df_plot_raw_rates)+i] <- df_plot_raw_rates[[i]]$Counts[3]/(df_plot_raw_rates[[i]]$Counts[1]+df_plot_raw_rates[[i]]$Counts[3])
# }
# 
# # get TP+FP and TP+FN
# sum_FP_TP <- df_plot[c(seq(1,60,3),seq(2,60,3)),"sum"]-rep(df_plot$Counts[seq(3,60,3)],2)
# sum_FN_TP <- df_plot[c(seq(3,60,3),seq(2,60,3)),"sum"]-rep(df_plot$Counts[seq(2,60,3)],2)
# 
# # extract TP and FPs
# df_plot <- df_plot[c(seq(1,60,3),seq(2,60,3)),]
# df_plot$prop <- c(prop1,1-prop1)
# df_plot$prop_min <- c(rep(0,20),prop1) 
# df_plot$prop_max <- c(prop1,rep(1,20)); rm(prop1)
# df_plot$sum <- sum_FP_TP ; rm(sum_FP_TP)
# df_plot$pie_area <- pi*(df_plot$sum/max(df_plot$sum))
# 
# # plot summary pies
# ggplot(df_plot, aes(fill=label, color=label))+ theme_classic() + facet_wrap(sample_id1~sample_id2,nrow = 5, ncol=5, drop = F) +
#   geom_rect(aes(xmin=prop_min,xmax=prop_max,ymin=0,ymax=pie_area),alpha=0.25,stat = "identity",na.rm = T) +   coord_polar("x") +
#   theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),legend.text=element_text(size=8), axis.ticks = element_blank(), strip.background = element_blank(),
#         strip.text.x = element_blank()) +
#   scale_fill_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                           "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   scale_color_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                            "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   ggsave(file = paste0(path, "Methylation/mb_mb_vs_mt_mt_DMPs/summary_figure_TP_FP.pdf"),device = "pdf", width=5, height = 5.5, units = "in")
# 
# df_plot <- df_plot[c(seq(3,60,3),seq(2,60,3)),]
# df_plot$prop <- c(prop2,1-prop2)
# df_plot$prop_min <- c(rep(0,20),prop2)
# df_plot$prop_max <- c(prop2,rep(1,20))
# df_plot$sum <- sum_FN_TP ; rm(sum_FN_TP)
# df_plot$pie_area <- pi*(df_plot$sum/max(df_plot$sum))
# 
# # plot summary pies
# ggplot(df_plot, aes(fill=label, color=label))+ theme_classic() + facet_wrap(sample_id1~sample_id2,nrow = 5, ncol=5, drop = F) +
#   geom_rect(aes(xmin=prop_min,xmax=prop_max,ymin=0,ymax=pie_area),alpha=0.25,stat = "identity",na.rm = T) +   coord_polar("x") +
#   theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),legend.text=element_text(size=8), axis.ticks = element_blank(), strip.background = element_blank(),
#         strip.text.x = element_blank()) +
#   scale_fill_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                           "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   scale_color_manual(name = "", values = c("true positives"="lightblue","false negatives"="lightsalmon","false positives"="darkred",
#                                            "true positives raw"="chartreuse3","false negatives raw"="mediumpurple4","false positives raw"="purple"),na.value="white") +
#   ggsave(file = paste0(path, "Methylation/mb_mb_vs_mt_mt_DMPs/summary_figure_TP_FN.pdf"),device = "pdf", width=5, height = 5.5, units = "in")
# 
# # Plot sensitivity and speciticity
# df_plot_rates = bind_rows(lapply(1:10, function(i,list) {
#   df <- list[[i]]
#   list[[i]]$TP <- df[df$label=="true positives","Counts"]; list[[i]]$FN <- df[df$label=="false negatives","Counts"]
#   list[[i]]$FP <- df[df$label=="false positives","Counts"]; list[[i]]$TN <- df[df$label=="true negatives","Counts"]
#   list[[i]]$sensitivity <- list[[i]]$TP/(list[[i]]$TP+list[[i]]$FN)
#   list[[i]]$accuracy <- list[[i]]$TP/(list[[i]]$TP+list[[i]]$FP)
#   # too many true negatives - specificity is not appropriate here. instead, we use the "accuracy" TP/TP+FP "accuracy" 
#   return(list[[i]]) }, list=df_plot_corr_rates))
# 
# df_plot_rates$labels <- unique(paste0(df_plot_corr$x[1:30],"_vs_",df_plot_corr$y[1:30]))
# df_plot_rates<- df_plot_rates[!duplicated(df_plot_rates[,3:6]),3:11]
# ggplot(data = df_plot_rates, aes(x = 1-sensitivity, y = 1-accuracy, label=labels)) + geom_point(shape=20,size=3) + geom_smooth(formula=y~exp(-x), method="lm") + theme_classic() + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) + 
#   geom_text(aes(x=1-sensitivity+0.005, label=labels),hjust=0, vjust=0, size = 3)
# 
# 
# # Plot sensitivity and speciticity
# df_plot_rates = bind_rows(lapply(1:10, function(i,list) {
#   df <- list[[i]]
#   list[[i]]$TP <- df[df$label=="true positives","Counts"]; list[[i]]$FN <- df[df$label=="false negatives","Counts"]
#   list[[i]]$FP <- df[df$label=="false positives","Counts"]; list[[i]]$TN <- df[df$label=="true negatives","Counts"]
#   list[[i]]$sensitivity <- list[[i]]$TP/(list[[i]]$TP+list[[i]]$FN)
#   list[[i]]$specificity <- list[[i]]$TN/(list[[i]]$TN+list[[i]]$FP) 
#   list[[i]]$accuracy <- list[[i]]$TP/(list[[i]]$TP+list[[i]]$FP)
#   list[[i]]$accuracy2 <- list[[i]]$TN/(list[[i]]$TN+list[[i]]$FN)
#   # too many true negatives - specificity is not appropriate here. instead, we use the "accuracy1" TP/TP+FP "accuracy1" 
#   return(list[[i]]) }, list=df_plot_raw_rates))
# 
# df_plot_rates$labels <- unique(paste0(df_plot_raw$x[1:30],"_vs_",df_plot_raw$y[1:30]))
# df_plot_rates<- df_plot_rates[!duplicated(df_plot_rates[,3:6]),3:11]
# ggplot(data = df_plot_rates, aes(x = sensitivity, y = accuracy, label=labels)) + geom_point(shape=20,size=3) + geom_smooth(formula=y~exp(-x), method="lm") + theme_classic() + coord_cartesian(xlim=c(0,1), ylim=c(0,1)) + 
#   geom_text(aes(x=sensitivity+0.005, label=labels),hjust=0, vjust=0, size = 3)

# Get all DMPs
DMPs <- DMPs_mb[,1:3]
tmp <- rowSums(DMPs_mt[,cols_corrected]==T)
DMPs$DMP_corrected_mt <- ifelse(tmp>1,T,F)

tmp <- rowSums(DMPs_mb[,cols_corrected]==T)
DMPs$DMP_corrected_mb <- ifelse(tmp>1,T,F)

tmp <- rowSums(DMPs_mt[,cols_raw]==T)
DMPs$DMP_raw_mt <- ifelse(tmp>1,T,F)

tmp <- rowSums(DMPs_mb[,cols_raw]==T)
DMPs$DMP_raw_mb <- ifelse(tmp>1,T,F)
rm(tmp)

df <- full_join(df, DMPs)
id <- c(grep("^m_n_R", colnames(df_DMPs)))
df$m_n <- rowMeans(df[,id], na.rm=T)

n = length(sample_ids); l = n-2 ; m = n+2
id <- c(grep("^m_t_corr_R", colnames(df)), which(colnames(df)=="m_n"))
df_DMPs <- df[!is.na(df$DMP_corrected_mt)&df$DMP_corrected_mt==T,id]
result_mt <- pvclust(df_DMPs, method.dist="cor", method.hclust="average", nboot=1000, r=seq(l/n,m/n,by=1/n))

id <- c(grep("^m_b_R", colnames(df)), which(colnames(df)=="m_n"))
df_DMPs <- df[!is.na(df$DMP_corrected_mb)&df$DMP_corrected_mb==T,id]
result_mb <- pvclust(df_DMPs, method.dist="cor", method.hclust="average", nboot=1000, r=seq(l/n,m/n,by=1/n))

dir.create(paste0(path,"/Clustering"), recursive = T)
pdf(file = paste0(path, "/Clustering/mt1_mt2_DMPs_clustering.pdf"))
plot(result_mt)
dev.off()
pdf(file = paste0(path, "/Clustering/mb1_mb2_DMPs_clustering.pdf"))
plot(result_mb)
dev.off()

# Executing the algorithm on curated data
perp = (nrow(t(df[df$DMP_corrected_mb==T,id]))/3)-1

id <- c(grep("^m_b_R", colnames(df)), which(colnames(df)=="m_n"))
tmp <- t(df[!is.na(df$DMP_corrected_mb)&df$DMP_corrected_mb==T,id])
umap_mb <- umap(tmp[,sample(1:ncol(tmp),5E4)])
col <- ifelse(rowSums(DMPs_mb[!is.na(df$DMP_corrected_mb)&df$DMP_corrected_mb==T,grepl("corrected_*R1",colnames(DMPs_mb))])>0,"red",
              ifelse(rowSums(DMPs_mb[!is.na(df$DMP_corrected_mb)&df$DMP_corrected_mb==T,grepl("corrected_R[1-4]_R5",colnames(DMPs_mb))])>0,"orange","white"))
ggplot(data.frame(umap_mb$layout))+geom_point(aes(x=X1,y=X2),shape=".",color=col) + theme_classic()

tsne <- Rtsne(t(tmp), dims = 2, perplexity=perp, verbose=TRUE, max_iter = 500)
tsne_mb <- data.frame(tsne$Y)
tsne_mb$label <- colnames(df[,id])

id <- c(grep("^m_t_corr_R", colnames(df)), which(colnames(df)=="m_n"))
tmp <- t(df[!is.na(df$DMP_corrected_mt)&df$DMP_corrected_mt==T,id])
tsne <- Rtsne(tmp, dims = 2, perplexity=perp, verbose=TRUE, max_iter = 500)
tsne_mt <- data.frame(tsne$Y)
tsne_mt$label <- colnames(df[,id])

tmp <- unlist(tmp)
tmp <- data.frame(m_t = as.numeric(tmp), labels = substr(names(tmp),1, str_count("m_x_Ri")))
umap_mb <- umap(tmp[,1])
rm(tsne, id) 

# keeping original data
df_tsne_mt_original=tsne_mt
df_tsne_mb_original=tsne_mb

# 3 cluster (normal, R1,R2,R3 and R4,R5)
nk=3

# Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
fit_cluster_hierarchical=hclust(dist(scale(t(df[!is.na(df$DMP_corrected_mt)&df$DMP_corrected_mt==T,id]))))
fit_cluster_hierarchical_bulk=hclust(dist(scale(t(df[!is.na(df$DMP_corrected_mb)&df$DMP_corrected_mb==T,id]))))

# setting nk clusters as output
df_tsne_mt_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=nk))
df_tsne_mb_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical_bulk, k=nk))  

plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="X1", y="X2", color=var_cluster)) +
    geom_point(size=3) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("tsne1") + ylab("tsne2") +
    ggtitle("") + theme_light(base_size=20) + theme_classic() +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.box = "horizontal",
          axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    geom_text(aes(label=ifelse(grepl("_n",label)==T,"normal",
                               substr(label,str_count("m_t_corr_")+1,str_count(label)))),hjust=0,size=3,col="black",
              nudge_x=(max(data$X1)-min(data$X1))/50,nudge_y=(max(data$X2)-min(data$X2))/50) +
    scale_colour_brewer(palette = palette)+
    coord_cartesian(xlim=c(min(data$X1)-abs(0.1*min(data$X1)),max(data$X1)+0.5*max(data$X1)))
}

# Plot side by side with gridExtra 
plot_k=plot_cluster(df_tsne_mt_original, "cl_hierarchical", "Accent")
plot_h=plot_cluster(df_tsne_mb_original, "cl_hierarchical", "Set1")
grid.arrange(plot_k, plot_h,  ncol=2)

grid.newpage()
df_summary <- data.frame("tx.ty" = substr(colnames(DMPs_mt)[cols_corrected],11,15), 
                         "common_SNV_clones" = c(4, 4, 1, 1, 5, 1, 1, 1, 1, 1.9),
                         num_DMPs_mtx_mty = as.numeric(colSums(na.omit(DMPs_mt[,cols_raw]==T))),
                         num_DMPs_mbx_mby = as.numeric(colSums(na.omit(DMPs_mb[,cols_raw]==T))),
                         num_DMPs_mtx_mty_corr = as.numeric(colSums(na.omit(DMPs_mt[,cols_corrected]==T))),
                         num_DMPs_mbx_mby_corr = as.numeric(colSums(na.omit(DMPs_mb[,cols_corrected]==T))))
plot(x=exp(df_summary$common_SNV_clones), y=df_summary$num_DMPs_mtx_mty)
points(x=exp(df_summary$common_SNV_clones), y=df_summary$num_DMPs_mbx_mby, col="red")

cols <- grep("_R1_R2",colnames(DMPs_mt))
DMPs <- na.omit(cbind.data.frame(DMPs_mt[,cols], DMPs_mb[,cols]))
grid.newpage()
venn.plot <- draw.quad.venn(
  area1 = sum(DMPs[,1]),
  area2 = sum(DMPs[,2]),
  area3 = sum(DMPs[,3]),
  area4 = sum(DMPs[,4]),
  n12 =  sum(DMPs[,1]==T&DMPs[,2]==T),
  n13 =  sum(DMPs[,1]==T&DMPs[,3]==T),
  n14 =  sum(DMPs[,1]==T&DMPs[,4]==T),
  n23 =  sum(DMPs[,2]==T&DMPs[,3]==T),
  n24 =  sum(DMPs[,2]==T&DMPs[,4]==T),
  n34 =  sum(DMPs[,3]==T&DMPs[,4]==T),
  n123 = sum(DMPs[,1]==T&DMPs[,2]==T&DMPs[,3]==T),
  n124 = sum(DMPs[,1]==T&DMPs[,2]==T&DMPs[,4]==T),
  n134 = sum(DMPs[,1]==T&DMPs[,3]==T&DMPs[,4]==T),
  n234 = sum(DMPs[,2]==T&DMPs[,3]==T&DMPs[,4]==T),
  n1234 =sum(DMPs[,1]==T&DMPs[,2]==T&DMPs[,3]==T&DMPs[,4]==T),
  category = c("mt1/mt2\nraw DMPs", "mt1/mt2\ncorrected DMPs", "mb1/mb2\nraw DMPs", "mb1/mb2\ncorrected DMPs"),
  col = c("orange", "red", "cyan", "blue"),
  lty = "solid",
  cex = 1,
  fontfamily = "Liberation Sans",
  cat.cex = 1,
  cat.fontfamily = "Liberation Sans", 
  cat.fontface = "bold",
  cat.col = c("orange", "red", "cyan", "blue")
)


#################################################
#################################################
####                                         ####
####                                         ####
####              tumour-tumour              ####
####                real data                ####
####                  v2                     ####
####                                         ####
#################################################
#################################################

source("~/.rprofile")
patient_info <- as.data.table(read.delim(paste0("~/rdeocd/txMeth/",
      "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"), 
      header=TRUE, sep="\t"))
patient_info[, group_rho := ifelse(rho_reseg < 0.3, 1, 
                            ifelse(rho_reseg > 0.58, 4, 0))]
table(patient_info$group_rho)
patient_ids <- unique(patient_info$Patient)

dt <- read.fst("~/rdeocd/txMeth/CpG_wise_CAMDAC_results_noNAs.fst",
               columns= c("ID","patient","sample","chromosome","start","end",
                          "m_t_low", "m_t_high", "m_t","m_b_low", "m_b_high", "m_b", "CN"),
               as.data.table=TRUE)
dt <- dt[chromosome!="X"]
dt[, LTXID := paste(patient,sample, sep="_")]
dt[, group_rho := patient_info$group_rho[match(LTXID, patient_info$LTX_ID)]]
dt[, rho := patient_info$rho_reseg[match(LTXID, patient_info$LTX_ID)]]
dt <- dt[, .SD, 
         .SDcols=c("ID","LTXID","patient","chromosome","start","end",
                   "m_t_low", "m_t_high", "m_t","m_b_low", "m_b_high", "m_b", 
                   "CN", "group_rho", "rho")]

# keep patients with low or high purity only
ids <- patient_info$LTX_ID[patient_info$group_rho==1|patient_info$group_rho==4]
dt <- dt[LTXID%in%ids]
rm(ids)

# set static params
set.seed(1234)
flags = c("intra", "inter")
CNs <- seq(1,6,1)
CN1 <- NULL
grs <- c(1,4)
reseg=TRUE
effect_size=0.2

for(flag in flags){
  # create results object
  dt_list <- list()
  # for(CN1 in CNs){
  for(gr1 in grs){
    # for(CN2 in CNs[CNs>=CN1]){
      for(gr2 in grs){
        if(!is.null(CN1)){
        print(paste0("copy number ", CN1, " vs. ", CN2 ,
                     " and purity ", gr1, " vs. ", gr2))
        } else {
          print(paste0("Purity groups ", gr1, " vs. ", gr2))
        }
        # get patient ids in purity group
        ids <- patient_info$LTX_ID[patient_info$group_rho==gr1|patient_info$group_rho==gr2]
        dt_tmp <- dt[LTXID%in%ids]
        rm(ids)
        
        if(!is.null(CN1)){
          # get CpG with desired copy numbers and purity combination
          if(CN2!=6&CN1!=6){
            dt_tmp <- dt_tmp[(CN==CN1&group_rho==gr1)|(CN==CN2&group_rho==gr2), ]
            if(flag=="inter"){
            idxs <- dt_tmp[, .(s1=sum(CN==CN1&!duplicated(patient)),
                               s2=sum(CN==CN2&!duplicated(patient))), by=.(ID)]
            } else {
            idxs <- dt_tmp[, .(s1=sum(CN==CN1&(duplicated(patient)|duplicated(patient, fromLast=TRUE))),
                               s2=sum(CN==CN2&(duplicated(patient)|duplicated(patient, fromLast=TRUE)))), by=.(ID)]
            }
          }    
          if(CN2!=6&CN1==6){
            dt_tmp <- dt_tmp[(CN>=CN1&group_rho==gr1)|(CN==CN2&group_rho==gr2), ]
            if(flag=="inter"){
            idxs <- dt_tmp[, .(s1=sum(CN>=CN1&!duplicated(patient)),
                               s2=sum(CN==CN2&!duplicated(patient))), by=.(ID)]
            } else {
            idxs <- dt_tmp[, .(s1=sum(CN>=CN1&(duplicated(patient)|duplicated(patient, fromLast=TRUE))),
                               s2=sum(CN==CN2&(duplicated(patient)|duplicated(patient, fromLast=TRUE)))), by=.(ID)]
            }
          }    
          if(CN2==6&CN1!=6){
            dt_tmp <- dt_tmp[(CN==CN1&group_rho==gr1)|(CN>=CN2&group_rho==gr2), ]
            if(flag=="inter"){
            idxs <- dt_tmp[, .(s1=sum(CN==CN1&!duplicated(patient)),
                               s2=sum(CN>=CN2&!duplicated(patient))), by=.(ID)]
            } else {
            idxs <- dt_tmp[, .(s1=sum(CN==CN1&(duplicated(patient)|duplicated(patient, fromLast=TRUE))),
                               s2=sum(CN>=CN2&(duplicated(patient)|duplicated(patient, fromLast=TRUE)))), by=.(ID)]
            }
          }    
          if(CN2==6&CN1==6){
            dt_tmp <- dt_tmp[(CN>=CN1&group_rho==gr1)|(CN>=CN2&group_rho==gr2), ]
            if(flag=="inter"){
            idxs <- dt_tmp[, .(s1=sum(CN>=CN1&!duplicated(patient)),
                               s2=sum(CN>=CN2&!duplicated(patient))), by=.(ID)]
            } else {
            idxs <- dt_tmp[, .(s1=sum(CN>=CN1&(duplicated(patient)|duplicated(patient, fromLast=TRUE))),
                               s2=sum(CN>=CN2&(duplicated(patient)|duplicated(patient, fromLast=TRUE)))), by=.(ID)]
            }
          }
        } else {
          dt_tmp <- dt_tmp[(group_rho==gr1)|(group_rho==gr2), ]
          if(flag=="inter"){
            idxs <- dt_tmp[, .(s1=sum(!duplicated(patient))), by=.(ID)]
          } else {
            idxs <- dt_tmp[, .(s1=sum(duplicated(patient)|duplicated(patient, fromLast=TRUE))), by=.(ID)]
          }
        }
        
        # select CpGs with info from different patients
        if(is.null(CN1)){
          idxs <- idxs[s1>=2, (ID)]
        } else {
          if(CN2==CN1){
            idxs <- idxs[s1>=2, (ID)]
          } else {
            idxs <- idxs[s1>=1&s2>=1, (ID)]
          }
        }
        dt_tmp <- dt_tmp[ID%in%idxs,]
        dt_tmp[, ID2 := paste(ID, patient, sep="_")]
        if(flag=="intra"){ 
          dt_tmp[, ID3 := paste(ID, LTXID, sep="_")]
          }
        rm(idxs)
        
        # separate remaining CpGs by CN and rho
        if(!is.null(CN1)){
          if(CN2!=6&CN1!=6){
            dt1 <- dt_tmp[CN==CN1&group_rho==gr1,]
            dt2 <- dt_tmp[CN==CN2&group_rho==gr2,]
          } 
          if(CN2!=6&CN1==6){
            dt1 <- dt_tmp[CN>=CN1&group_rho==gr1,]
            dt2 <- dt_tmp[CN==CN2&group_rho==gr2,]
          }    
          if(CN2==6&CN1!=6){
            dt1 <- dt_tmp[CN==CN1&group_rho==gr1,]
            dt2 <- dt_tmp[CN>=CN2&group_rho==gr2,]
          }    
          if(CN2==6&CN1==6){
            dt1 <- dt_tmp[CN>=CN1&group_rho==gr1,]
            dt2 <- dt_tmp[CN>=CN2&group_rho==gr2,]
          }
        } else {
          dt1 <- dt_tmp[group_rho==gr1,]
          dt2 <- dt_tmp[group_rho==gr2,]
        }
        
        # get CpGs unique IDs
        if(flag=="inter"){ 
        dt1 <- dt1[!duplicated(ID),]
        dt1 <- dt1[!duplicated(ID2),] # pick one entry per patient max for each CpG in dt1
        dt2 <- dt2[!ID2%in%dt1$ID2,] # pick all CpGs not already selected for a given patient in dt1
        dt2 <- dt2[!duplicated(ID)]
        
        # match IDs
        dt2 <- dt2[ID%in%dt1$ID]
        dt1 <- dt1[ID%in%dt2$ID]
        
        } else {
        dt1 <- dt1[!duplicated(ID2),] # pick one entry per patient max for each CpG in dt1
        dt2 <- dt2[(!ID3 %in% dt1$ID3),] # pick all CpGs not in dt1 
        dt2 <- dt2[!duplicated(ID2)] # pick one entry per patient max for each CpG not in dt1
        
        # match IDs
        dt2 <- dt2[ID2%in%dt1$ID2,]
        dt1 <- dt1[ID2%in%dt2$ID2,]
        
        # assign patient + CpG coord to ID column
        dt1[, ID := ID2,]
        dt2[, ID := ID2,]
        }
        
        cols <- c("ID","m_t_low", "m_t_high", "m_t","m_b_low", "m_b_high", "m_b")
        dt1[, c("CN1", "rho1", "id1", "gr1") := list(CN, rho, LTXID, group_rho)]
        dt2[, c("CN2", "rho2", "id2", "gr2") := list(CN, rho, LTXID, group_rho)]
        dt1 <- dt1[, .SD, .SDcols=c(cols,"CN1", "rho1", "id1", "gr1")]
        dt2 <- dt2[, .SD, .SDcols=c(cols,"CN2", "rho2", "id2", "gr2")]
          
        cols1 <- paste(colnames(dt1)[grepl("m_", colnames(dt1))],1,sep="_")
        cols2 <- paste(colnames(dt2)[grepl("m_", colnames(dt2))],2,sep="_")
        colnames(dt1)[grepl("m_", colnames(dt1))] <- cols1 
        colnames(dt2)[grepl("m_", colnames(dt2))] <- cols2   
        rm(cols1, cols2)
        
        df <- Reduce(function(...) merge(..., by=c("ID"),all = T), list(dt1, dt2))
        dt_tmp <- data.table(df)
        
        # print stats
        cat(nrow(dt_tmp), " CpGs to sample from", "\n")
        
        # sample
        sample_size <- ifelse(!is.null(CN1), 1E4, 1E6)
        if(nrow(dt_tmp)>=1E3){
        if(nrow(dt_tmp)<1E4){
          dt_tmp <- dt_tmp[sample(1:nrow(dt_tmp), sample_size, replace=TRUE)]
        } else {
          dt_tmp <- dt_tmp[sample(1:nrow(dt_tmp), sample_size, replace=FALSE)]
        }
        rm(df)
        
        # get overlaps
        cols1 <- colnames(dt1)[grepl("m_t_[hl]", colnames(dt1))]
        cols2 <- colnames(dt2)[grepl("m_t_[hl]", colnames(dt2))]
        raw <- Overlap(dt_tmp[, .SD, .SDcols=cols1],dt_tmp[, .SD, .SDcols=cols2])==0
        diffs <- dt_tmp[, m_t_1-m_t_2] 
        corr <- abs(diffs) >= effect_size & raw == T
        rm(diffs)
        
        # store
        dt_tmp[, DMP_tt_raw := raw]
        dt_tmp[, DMP_tt_corr := corr]
        
        # get overlaps
        cols1 <- colnames(dt1)[grepl("m_b_[hl]", colnames(dt1))]
        cols2 <- colnames(dt2)[grepl("m_b_[hl]", colnames(dt2))]
        raw <- Overlap(dt_tmp[, .SD, .SDcols=cols1],dt_tmp[, .SD, .SDcols=cols2])==0
        diffs <- dt_tmp[, m_b_1-m_b_2] 
        corr <- abs(diffs) >= effect_size & raw == T
        rm(diffs)
        
        # store
        dt_tmp[, DMP_bb_raw := raw]
        dt_tmp[, DMP_bb_corr := corr]
        dt_list[[length(dt_list)+1]] <- dt_tmp
        rm(raw,corr,cols1,cols2,dt1,dt2)
        }
      # }
    }
  # }
}

    
 # get all DMPs loci DMPs between at least one pair in bulk or pure tumour-tumour DMP set
 dt_results <- rbindlist(dt_list)
 assign(paste0("dt_list_", flag), dt_list)
 cols_corrected <- grep(("_corr"),colnames(dt_results))
 rm(dt_list)

 # set vars
 venn_corr <- list(); venn_raw <- list()
 max_DMPs_corr <- 1E4 
 # dir_path <- paste0("~/rdeocd/txMeth_analysis/mb_mb_vs_mt_mt_DMPs/venn_plots/")
 # dir.create(dir_path)
 
 dt_results[ , labels := ifelse(gr1 == 1 & gr2 == 1, "low_vs_low",
                         ifelse(gr1 == 4 & gr2 == 4, "high_vs_high", "low_vs_high"))]
 
 labs <- unique(dt_results$labels)
 for(label in labs){
 # for(CN1 in CNs){
   # for(group1 in grs){
     # for(CN2 in CNs){
       # for(group2 in grs){
         tmp <- dt_results[labels == label,]
         cols <- colnames(tmp)[grepl("DMP_[tb][tb]_corr",colnames(tmp))] 
         
          # Get venn counts (corrected calls)
          pdf(paste0(out_dir,flag,"_",label,"_mt_mt_vs_mb_mb_DMPs_corrected_new.pdf"))
          print(plot(euler(tmp[,.SD, .SDcols=cols], shape = "ellipse"), quantities = TRUE))
          dev.off()
          Sys.sleep(1)
      # }
    # }
  # }
 }
}
save(dt_list_inter, dt_list_intra, file="~/rdeocd/txMeth/power_real_results_list.RData")
rm(list=ls(pattern="CN|gr"))

# rm(df.venn, df.vdc.ordered, tmp, cols, coordinates, x_coord, y_coord, d)
# 
# # merge dfs for plotting
# df_plot <-  bind_rows(venn_corr) 
# rm(tmp)
# 
# df_plot$gr1 <- factor(df_plot$gr1, levels=seq(1,4,1), ordered = TRUE)
# df_plot$gr2 <- factor(df_plot$gr2, levels=seq(1,4,1), ordered = TRUE)
# 
# # Plot summary venn diagram
# ggplot(data=df_plot[is.finite(df_plot$x)&grepl("raw", df_plot$labels),]) + coord_fixed() + theme_void() + 
#   facet_wrap(gr1~gr2, ncol = length(sample_ids2), nrow=length(sample_ids1), drop = F) +
#   geom_circle(aes(x0 = x, y0 = y, r = rad, color = labels, fill = labels), alpha = .25, size=.5) +  
#   theme(legend.position = 'right',axis.text = element_blank(), axis.title = element_blank(), 
#         axis.line = element_blank(),legend.text=element_text(size=8), 
#         axis.ticks = element_blank(), strip.background = element_blank(),
#         #strip.text.x = element_blank()
#   ) +
#   scale_colour_manual(name="legend", values = c("mb1/mb2\n DMPs"="lightgrey","mt1/mt2\n DMPs"="dodgerblue",
#                                                 "mb1/mb2\n raw DMPs"="grey20","mt1/mt2\n raw DMPs"="darkblue")) + 
#   scale_fill_manual(name="legend",values = c("mb1/mb2\n DMPs"="lightgrey","mb1/mb2\n raw DMPs"="grey20",
#                                              "mt1/mt2\n DMPs"="dodgerblue","mt1/mt2\n raw DMPs"="darkblue")) + 
#   ggsave(file = paste0(out_dir, "/summary_figure_venn_raw.pdf"),device = "pdf", width=5, height = 6.5, units = "in") 
# 
# # Plot summary venn diagram
# ggplot(data=df_plot[is.finite(df_plot$x)&!grepl("raw", df_plot$labels),]) + coord_fixed() + theme_void() + 
#   facet_wrap(sample_id1~sample_id2, ncol = length(sample_ids2), nrow=length(sample_ids1), drop = F) +
#   geom_circle(aes(x0 = x, y0 = y, r = rad, color = labels, fill = labels), alpha = .25, size=.5) +  
#   theme(legend.position = 'right',axis.text = element_blank(), axis.title = element_blank(), 
#         axis.line = element_blank(),legend.text=element_text(size=8), 
#         axis.ticks = element_blank(), strip.background = element_blank(),
#         #strip.text.x = element_blank()
#   ) +
#   scale_colour_manual(name="legend", values = c("mb1/mb2\n DMPs"="lightgrey","mt1/mt2\n DMPs"="dodgerblue",
#                                                 "mb1/mb2\n raw DMPs"="grey20","mt1/mt2\n raw DMPs"="darkblue")) + 
#   scale_fill_manual(name="legend",values = c("mb1/mb2\n DMPs"="lightgrey","mb1/mb2\n raw DMPs"="grey20",
#                                              "mt1/mt2\n DMPs"="dodgerblue","mt1/mt2\n raw DMPs"="darkblue")) + 
#   ggsave(file = paste0(out_dir, "/summary_figure_venn_corrected.pdf"),device = "pdf", width=5, height = 6.5, units = "in") 


