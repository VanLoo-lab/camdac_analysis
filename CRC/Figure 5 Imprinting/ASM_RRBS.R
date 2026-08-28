setwd("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting")

source("./src/hetSNP_phase.R")
x <- c("Biobase","GenomicRanges","Rsamtools","data.table","stringr","readr","dplyr","matrixStats","ggplot2","gridExtra", "parallel", "doParallel")
invisible(lapply(x, function(y) suppressWarnings(suppressMessages(library(y, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))))
rm(x)

# load the sample information and run CAMDAC_ASM for RRBS -------
CN_good <- readLines("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/Good_samples_59.txt")
metadata <- read.csv("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/CAMDAC_metadata_MSS.csv")
meta_rrbs <- read.csv("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/CN_103/rrbs_cn_summary.csv")
metadata <- inner_join(metadata, meta_rrbs, by=c("PT_ID"="sample_id"))
metadata <- metadata %>% filter(PT_ID %in% CN_good)
loci_file = "/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/loci.hg38.all.RData"
n_cores=16
mq=0
build = "hg38"
path_to_CAMDAC <- "~/R/x86_64-pc-linux-gnu-library/4.1/CAMDAC/"

# load reference file
load(loci_file)

for(jobid in 1:length(CN_good)){
    tryCatch({
    print(paste0("Working on ", jobid, " of the 59 patients..."))

    PT_ID <- metadata$PT_ID[jobid]
    sex <- metadata$Sex[jobid]
    sample_ids = c("N", "T")
    path = file.path("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS",PT_ID)

    for(sample_id in sample_ids){
    # Load methylation data and get bam_file path 
    if(sample_id=="N"){
        print("Working on normal allele counting...")
        bam_file <- metadata$bam_normal[jobid]
        load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Methylation/N1/dt_normal_m.RData"))
        gr <- makeGRangesFromDataFrame(dt_normal_m)
        rm(dt_normal_m)
        
        #Get Het SNPs
        load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Copy_number/N1/",PT_ID,".N1.SNPs.RData"))
        df <- dt[dt$BAF>=0.15&dt$BAF<=0.85,]
        rm(dt)
    }

    if(sample_id=="T"){
        print("Working on tumor allele counting...")
        bam_file <- metadata$bam_tumor[jobid]
        load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Methylation/T1/CAMDAC_results_per_CpG.RData"))
        gr <- makeGRangesFromDataFrame(CAMDAC_results_per_CpG)
        seqlevelsStyle(gr) <- "UCSC"
        rm(CAMDAC_results_per_CpG)
        
        #Get Het SNPs
        load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Copy_number/T1/",PT_ID,".T1.SNPs.RData"))
        df <- dt[dt$BAF_n>=0.15&dt$BAF_n<=0.85,]
        rm(dt) 
    }

    # Extract het SNPs positions
    df_POS <- data.table(df[,c("chrom","POS","ref","alt")])

    # Make het SNPs GRanges
    gr_SNPs <- GRanges(seqnames=df$chrom, IRanges(start=df$POS, end=df$POS))
    seqlevelsStyle(gr_SNPs) <- "UCSC"
    rm(df)

    # Get CpGs which could be within 1 read length to nearead het SNP
    distances <- data.frame(mcols(distanceToNearest(gr, gr_SNPs)))[,1]
    gr <- gr[distances<=120]
    rm(distances)

    # Set options and working directory
    suppressWarnings(dir.create(file.path(paste0(path, "/", sample_id)), recursive = TRUE))
    path_output=path
    setwd(path_output)

    # Set build variables
    if(build%in%c("hg19","hg38")){UCSC=TRUE}
    if(build%in%c("GRCH37","GRCH38")){UCSC=FALSE}
    if(build=="GRCH37"){build="hg19"} # set build to to assembly version disregarging USCS/Ensembl
    if(build=="GRCH38"){build="hg38"}
    cat(paste("ScanBam pileup with build ", build, sep = " "), "\n", sep = "")

    # Only keep het SNPs and nearby CpGs
    loci_subset <- all_loci[(all_loci%in%gr)|(all_loci%over%gr_SNPs)]
    names(loci_subset) <- 1:length(loci_subset)
        
    # get_reads in parrallel
    if(n_cores>1){
    print(paste0("Multiple cores detected: using ", n_cores, " cores in parallel..."))
    
    # Set the cluster
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Pass libPaths to each worker on the cluster if incorrect libPaths in Rprofile.
    if(file.exists("~/.rprofile")|file.exists("~/.Rprofile")){
        x <- .libPaths()
        suppressWarnings(try(source("~/.rprofile"), silent=TRUE))
        y <- .libPaths()
        if(sum(!x%in%y)>0){
        clusterCall(cl, function(x) .libPaths(x), .libPaths())
        }
        rm(x, y)
    }
    
    # run get reads function
    df_phased <- NULL
    df_phased <- foreach(j = 1:n_cores, .combine='rbind.data.frame', 
                            .packages=c("magrittr","Rsamtools", "GenomicRanges", "S4Vectors", "IRanges", "dplyr", "data.table"),
                            .multicombine = TRUE, .init=list(list(), list())) %dopar% {
                            df_Jth_core <- get_reads(i=i,j=j,n_cores=n_cores,bam_file=bam_file,loci_subset=loci_subset,df_POS=df_POS,test=FALSE, mq=0)
                            }
    stopCluster(cl)
    df_phased <- data.table(df_phased)
    }
    
    if(n_cores==1){df_merged <- get_reads(i=i,j=NULL,n_cores=NULL,bam_file=bam_file,loci_subset=loci_subset,df_POS=df_POS,test=FALSE, mq=0)}
    invisible(gc())
    
    # subset to simple genotypes
    df_phased <- df_phased[genotype%in%c("ref", "alt")]
    
    # in case of duplicate positions due to SNP, keep pos with methylation rate
    df_phased[, l := length(genotype), by=.(chrom,start,end)]
    df_phased[, l2 := length(unique(genotype)), by=.(chrom,start,end)]
    df_phased[, l3 := length(unique(POS_genotype)), by=.(chrom,start,end)]
    
    # filter df
    df_phased <- df_phased[l==2&l2==2&l3==1,]
    
    save(df_phased, file=paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/",PT_ID,"/",
                                            sample_id,"/allele_specific_allele_counts_",sample_id,".RData"))
    print("Allele counting object saved.")

    if(sample_id == "N"){ 
        # get allele specific methylation
        print("Generating ASM for normal...")
        df_phased_methylation <- df_phased[order(chrom,start,end,width,POS,genotype), 
                                    .(genotype_1="ref",genotype_2="alt", num_genotypes=length(genotype),
                                    POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                    total_depth_n_1 = total_depth[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                    total_depth_n_2 = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                    POS_genotype_BAF = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_depth),
                                    total_counts_n_1 = total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                    total_counts_n_2 = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                    CpG_BAF_n = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_counts),
                                    M_n_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                    UM_n_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                    total_counts_m_n_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                    M_n_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                    UM_n_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                    total_counts_m_n_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]), 
                                    keyby=c("CHR","chrom","start","end","width","POS","ref","alt")]  
        rm(df_phased)
        
        # trim loci with insufficient counts in normal 
        df_phased_methylation <- df_phased_methylation[(total_counts_m_n_1 >= 5 & total_counts_m_n_2 >= 5) | 
                                                        (!is.na(POS) &paste(chrom, POS) %in% paste(chrom, POS_genotype)), ]
        
        # get potentially imprinted loci
        df_phased_methylation[, m_n_1 := as.numeric(M_n_1)/as.numeric(total_counts_m_n_1)] 
        df_phased_methylation[, m_n_2 := as.numeric(M_n_2)/as.numeric(total_counts_m_n_2)]
        df_phased_methylation[, m_diff_n12 := abs(m_n_1-m_n_2)]
        
        # Calculate HDI for the methylation rate
        n <- nrow(df_phased_methylation)
        M=df_phased_methylation$M_n_1
        UM=df_phased_methylation$UM_n_1
        vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
        vec[,1:2]<- do.call(rbind,mclapply(1:n, function(i,M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M[i]+1,shape2=UM[i]+1),mc.cores=n_cores,M=M,UM=UM))
        df_phased_methylation[, as.character(c("m_n_1_low", "m_n_1_high")) := as.list(vec)]
        
        # Calculate HDI for the methylation rate
        M=df_phased_methylation$M_n_2
        UM=df_phased_methylation$UM_n_2
        vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
        vec[,1:2]<- do.call(rbind,mclapply(1:n, function(i,M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M[i]+1,shape2=UM[i]+1),mc.cores=n_cores,M=M,UM=UM))
        df_phased_methylation[, as.character(c("m_n_2_low", "m_n_2_high")) := as.list(vec)]
        
        # check overlaps
        df_phased_methylation[, is_DM_n12 :=  ((m_n_1_low>m_n_2_high)|(m_n_2_low>m_n_1_high)) & m_diff_n12 > 0.5]
        # df_phased_methylation[,table( is_DM_n12, round2(m_diff_n12,digits=1))]
        
        # get bulk methylation rates
        df_phased_methylation[, M_n := as.numeric(M_n_1)+as.numeric(M_n_2)] 
        df_phased_methylation[, UM_n := as.numeric(UM_n_1)+as.numeric(UM_n_2)]
        df_phased_methylation[, total_counts_m_n := M_n + UM_n]
        df_phased_methylation[, m_n := M_n/total_counts_m_n]
        
        # Calculate HDI for the methylation rate
        M=df_phased_methylation$M_n
        UM=df_phased_methylation$UM_n
        vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
        vec[,1:2]<- do.call(rbind,mclapply(1:n, function(i,M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M[i]+1,shape2=UM[i]+1),mc.cores=n_cores,M=M,UM=UM))
        df_phased_methylation[, as.character(c("m_n_low", "m_n_high")) := as.list(vec)]
        rm(vec, n, UM, M)
   
        save(df_phased_methylation, file=paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/",PT_ID,"/N/allele_specific_methylation_rates_N.RData"))
        print("ASM_N object saved.")
    }
        
    if(sample_id == "T"){
        # load normal, bulk and pure tumour methylation rates (no phasing)
        print("Generating ASM for tumor...")
        load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Methylation/T1/CAMDAC_results_per_CpG.RData"))
         
        # get columns of interest
        df_mt <- data.table(CAMDAC_results_per_CpG)
        dt1 <- df_mt[,.(chrom, start, end, nA, nB, label, "CpG_BAF"=BAF)]
        dt2 <- df_mt[,.(chrom, start, end, m_b, cov_b, m_b_low, m_b_high,
                                    m_t_corr, cov_t, m_t_low, m_t_high, DMP_t)]
        
        # annotate copy number and BAF on phased methylation rates  
        df_tumour <- merge(df_phased, dt1, by = c("chrom", "start", "end"), all=FALSE)
        df_tumour <- df_tumour[!is.na(nA)&!is.na(total_counts_m),]
        rm(dt1, df_phased, df_mt)
        
        # load matched normal
        fls = paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/",PT_ID,"/N/allele_specific_methylation_rates_N.RData")
        load(fls)
        
        # remove loci with no coverage in the normal
        df_normal <- df_phased_methylation[!is.na(total_counts_m_n), ]
        rm(df_phased_methylation, fls)
            
        # get allele-specific methylation
        df_tumour[, CN := nA+nB]
        df_phased_methylation <- df_tumour[order(chrom,start,end,width,POS,genotype), 
                                            .(genotype_1="ref", genotype_2="alt", num_genotypes=length(genotype),
                                            POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                            total_depth_1 = total_depth[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                            total_depth_2 = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                            POS_genotype_BAF_t = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_depth),
                                            total_counts_1 = total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                            total_counts_2 = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                            CpG_BAF_t = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_counts),
                                            M_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                            UM_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                            total_counts_m_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                            M_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                            UM_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                            total_counts_m_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]),
                                            keyby=c("CHR","chrom","start","end","width","POS","ref","alt", "nA", "nB", "CN")]  
        
        # label and order cols
        df_phased_methylation[, label := ifelse((POS==start & ref=="C")|(POS==end & ref=="G"), "CG-destroying","CG-forming")]
        df_phased_methylation <- 
        df_phased_methylation[, .SD, .SDcols=c(1:which(colnames(df_phased_methylation)=="CpG_BAF_t"), which(colnames(df_phased_methylation)=="label"),
                                                (which(colnames(df_phased_methylation)=="CpG_BAF_t")+1):which(colnames(df_phased_methylation)=="total_counts_m_2"))]
        
        # join tumour and normal info
        df_phased_methylation <- merge(df_phased_methylation, df_normal, 
                                        by = c("CHR", "chrom", "start", "end", "width","POS", "ref", "alt", 
                                        "POS_genotype", "genotype_1", "genotype_2", "num_genotypes"), all=FALSE)
        rm(df_normal)
        
        # get phased methylation rates
        df_phased_methylation[, m_1 := as.numeric(M_1)/as.numeric(total_counts_m_1)] 
        df_phased_methylation[, m_2 := as.numeric(M_2)/as.numeric(total_counts_m_2)]
        df_phased_methylation[ ,m_diff_b12 := abs(m_1-m_2)]
        
        # Calculate HDI for the methylation rate
        n <- nrow(df_phased_methylation)
        M=df_phased_methylation$M_1
        UM=df_phased_methylation$UM_1
        vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
        vec[,1:2]<- do.call(rbind,mclapply(1:n, function(i,M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M[i]+1,shape2=UM[i]+1),mc.cores=n_cores,M=M,UM=UM))
        df_phased_methylation[, as.character(c("m_1_low", "m_1_high")) := as.list(vec)]
        
        # Calculate HDI for the methylation rate
        M=df_phased_methylation$M_2
        UM=df_phased_methylation$UM_2
        vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
        vec[,1:2]<- do.call(rbind,mclapply(1:n, function(i,M,UM) HDIofICDF(ICDFname=qbeta,credMass=.99,shape1=M[i]+1,shape2=UM[i]+1),mc.cores=n_cores,M=M,UM=UM))
        df_phased_methylation[, as.character(c("m_2_low", "m_2_high")) := as.list(vec)]
        
        # check overlaps
        df_phased_methylation[, is_DM_b12 :=  ((m_1_low>m_2_high)|(m_2_low>m_1_high)) & m_diff_b12 > 0.5]
        
        # extract purity
        rho <- metadata$purity[jobid]
        rho <- as.numeric(unname(unlist(rho)))
        
        # # load phased haplotypes
        # fls = paste0(parent_dir,PT_ID,"/Copy_number/",PT_ID,"_haplotypes")
        # fls = paste0(fls,".RData")
        # load(fls)
        
        # # set for overlap 
        # dt_phasing <- data.table(dt_phasing)
        # dt_phasing[, c("POS_genotype") := list(as.character(POS))]
        # dt_phasing <- dt_phasing[, .SD, .SDcols=c("chrom", "POS_genotype", paste0("phased_",sample_id), "phasing", "re_seg")]
        # colnames(dt_phasing)[colnames(dt_phasing)==paste0("phased_",sample_id)] <- "phased_seg_BAF"
        
        # # overlap phased methylation rates with haplotypes
        # df_phased_methylation[, POS_genotype := tstrsplit(POS_genotype, ",", keep=1)]
        # tmp <- merge(df_phased_methylation, dt_phasing, by = c("chrom", "POS_genotype"), all.x=TRUE)
        # tmp[, phasing_tmp := ifelse(is.na(phasing), "unknown", phasing)]
        # tmp[, phasing := phasing_tmp]
        # tmp[, phasing_tmp := NULL]
        # rm(dt_phasing)
        
        ### this is added by lesky: if there is no haplotype object, just add a column "phasing" to be "unknown", test the remaining code
        tmp <- df_phased_methylation
        tmp$phasing <- "unknown"
        ### end of added lines

        # get purified tumour phased methylation rates
        df_phased_methylation <- get_phased_pure_tumour_methylome(tmp, p=rho, sex = sex)
        rm(tmp)
        
        #tmp <- get_DMPs(path=NULL,df=df_phased_methylation,n_cores=1)
        df_phased_methylation[, m_t_1_raw := m_t_1] 
        df_phased_methylation[, m_t_2_raw := m_t_2]
        df_phased_methylation[, m_t_1 := ifelse(m_t_1_raw<0, 0, 
                                                ifelse(m_t_1_raw>1, 1, m_t_1_raw))] 
        df_phased_methylation[, m_t_2 := ifelse(m_t_2_raw<0, 0, 
                                                ifelse(m_t_2_raw>1, 1, m_t_2_raw))]
            
        # get allelic bias
        df_phased_methylation[,m_diff_t12 := abs(m_t_1-m_t_2)]

        # get methylation difference
        df_phased_methylation[,m_diff_tn1 := abs(m_n_1-m_t_1)]
        df_phased_methylation[,m_diff_tn2 := abs(m_n_2-m_t_2)]
        df_phased_methylation[,m_diff_bn1 := abs(m_n_1-m_1)]
        df_phased_methylation[,m_diff_bn2 := abs(m_n_2-m_2)]
            
        # Compute the tumour methylation rate HDI 
        n <- nrow(df_phased_methylation)
        df_phased_methylation[, as.character(c("m_t_1_low", "m_t_1_high")) := 
                                as.list(HDIofMCMC(M_b=M_1,UM_b=UM_1,M_n=M_n_1,UM_n=UM_n_1,p=rho,CN=n1,CN_n=n1_n,0.99)), by=.(1:n)]
            
        # cap allele 1 HDI values
        df_phased_methylation[, m_t_1_low := ifelse(m_t_1_low < 0, 0, ifelse(m_t_1_low > 1 ,1, m_t_1_low))]
        df_phased_methylation[, m_t_1_high := ifelse(m_t_1_high < 0, 0, ifelse(m_t_1_high > 1 ,1, m_t_1_high))]
            
        # Compute the tumour methylation rate HDI 
        ## Allele 2
        df_phased_methylation[, as.character(c("m_t_2_low", "m_t_2_high")) := 
                                as.list(HDIofMCMC(M_b=M_2,UM_b=UM_2,M_n=M_n_2,UM_n=UM_n_2,p=rho,CN=n2,CN_n=n2_n,0.99)), by=.(1:n)]
        rm(n)
            
        # cap allele 2 HDI values
        df_phased_methylation[, m_t_2_low := ifelse(m_t_2_low < 0, 0, ifelse(m_t_2_low > 1 ,1, m_t_2_low))]
        df_phased_methylation[, m_t_2_high := ifelse(m_t_2_high < 0, 0, ifelse(m_t_2_high > 1 ,1, m_t_2_high))]
        
        # # Compute exact overlap
        # df_phased_methylation[, ov_t12 := do.call(
        #   rbind,mclapply(1:nrow(df_phased_methylation), 
        #                  function(i,M_b_1,UM_b_1,M_n_1,UM_n_1,CN_1,CN_n_1,
        #                           M_b_2,UM_b_2,M_n_2,UM_n_2,CN_2,CN_n_2,p) {
        #     MCMC_overlap_pure(M_b_1[i],UM_b_1[i],M_n_1[i],UM_n_1[i],CN_1[i],CN_n_1[i],
        #                       M_b_2[i],UM_b_2[i],M_n_2[i],UM_n_2[i],CN_2[i],CN_n_2[i],p)
        #   },mc.cores=n_cores,
        #   M_b_1=M_1,UM_b_1=UM_1,M_n_1=M_n_1,UM_n_1=UM_n_1,CN_1=n1,CN_n_1=n1_n,
        #   M_b_2=M_2,UM_b_2=UM_2,M_n_2=M_n_2,UM_n_2=UM_n_2,CN_2=n2,CN_n_2=n2_n,
        #   p=rho))]
            
        # check overlaps
        df_phased_methylation[, is_DM_t12 :=  ((m_t_1_low>m_t_2_high)|(m_t_2_low>m_t_1_high)) & m_diff_t12 > 0.5]
        #df_phased_methylation[!is.na(is_DM_t12)&ov_t12<0.21&m_diff_t12>0.5, table(is_DM_t12==TRUE,round2(ov_t12,digits=2))] 
            
        # Checkpoint
        print("Purified tumour methylation rate HDI added for each CpG locus per parental alleles")
            
        # # annotate imprinted gene loci
        # setkeyv(df_phased_methylation, cols=c("chrom", "start", "end"))
        # df_phased_methylation <- foverlaps(df_phased_methylation, ov)
        # df_phased_methylation <- df_phased_methylation[!is.na(CHR)&!is.na(Gene)]
                    
        # # Create file
        # f_nm <- paste0("~/rdeocd/txMeth/",PT_ID,"/Allele_specific_analysis/",
        #                sample_id,"/allele_specific_methylation_rates_",sample_id,".fst")
        # write_fst(df_phased_methylation,  f_nm)
            
        # join phased and bulk tumour info
        df_phased_methylation <- merge(df_phased_methylation, dt2, by = c("chrom", "start", "end"), all=FALSE)
        
        # order colnames
        df_phased_methylation <- 
        df_phased_methylation[, .SD, .SDcols=c(1:which(colnames(df_phased_methylation)=="label"),
                                                which(colnames(df_phased_methylation)=="total_counts_n_1"):which(colnames(df_phased_methylation)=="CpG_BAF_n"),
                                                which(colnames(df_phased_methylation)=="phasing"),which(colnames(df_phased_methylation)=="phased_seg_BAF"),
                                                which(colnames(df_phased_methylation)=="n1"):which(colnames(df_phased_methylation)=="n2_n"),
                                                which(colnames(df_phased_methylation)=="M_1"):which(colnames(df_phased_methylation)=="total_counts_m_2"),
                                                which(colnames(df_phased_methylation)=="m_1"):which(colnames(df_phased_methylation)=="is_DM_b12"),
                                                which(colnames(df_phased_methylation)=="m_b"):which(colnames(df_phased_methylation)=="m_b_high"),
                                                which(colnames(df_phased_methylation)=="M_n_1"):which(colnames(df_phased_methylation)=="m_n_high"), 
                                                which(colnames(df_phased_methylation)=="m_t_1"):which(colnames(df_phased_methylation)=="is_DM_t12"), 
                                                which(colnames(df_phased_methylation)=="m_t_corr"):which(colnames(df_phased_methylation)=="DMP_t"))]
        #save
        save(df_phased_methylation, file=paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/", PT_ID,
                                                "/T/allele_specific_methylation_rates_T.RData"))
        print("ASM_T object saved.")
    }
    }
   })
}  

rm(all_loci, df_phased_methylation, UCSC, sex, sample_id, sample_ids,
   PT_ID, mq, bam_file, gr, df_POS, path)




# ASM version 2 for potential imprinted genes -------
# Lili's code rename previous version to v2 and keep the output from the following as no "_v2" in file name
# my updated code keep the output from the following as "v1" and the original as "" in file name
for(jobid in 1:length(CN_good)){
    tryCatch({
    print(paste0("Working on ", jobid, " of the 59 patients for imprinting...")) 
    
    PT_ID <- metadata$PT_ID[jobid]
    sex <- metadata$Sex[jobid]
    sample_ids = c("T")
    path = file.path("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS",PT_ID)

    for(sample_id in sample_ids){ 
    # Load phased methylation information
    print(paste(PT_ID,sample_id,sep = "-"))
    fls <- paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/",PT_ID,"/",
                                            sample_id,"/allele_specific_allele_counts_",sample_id,".RData")
    load(fls)
    
    if(sample_id != "N"){
      # load normal, bulk and pure tumour methylation rates (no phasing)
      load(paste0("/rsrch9/scratch/genetics/vanloolab/Lesky/CRC_RRBS/", metadata$folder[jobid], "/bam_final/", PT_ID,"/Methylation/T1/CAMDAC_results_per_CpG.RData"))
      
      # get columns of interest
      df_mt <- data.table(CAMDAC_results_per_CpG)
      dt1 <- df_mt[,.(chrom, start, end, nA, nB)]
      dt2 <- df_mt[,.(chrom, start, end, m_b, cov_b, m_b_low, m_b_high,
                      m_t_corr, cov_t, m_t_low, m_t_high, DMP_t)]
      
      # annotate copy number and BAF on phased methylation rates  
      df_tumour <- merge(df_phased, dt1, by = c("chrom", "start", "end"), all=FALSE)
      df_tumour <- df_tumour[!is.na(nA)&!is.na(total_counts_m),]
      rm(dt1, df_phased, df_mt)
      
      # load matched normal
      fls = paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/",PT_ID,"/N/allele_specific_methylation_rates_N.RData")
      load(fls)
      
      # remove loci with no coverage in the normal
      df_normal <- df_phased_methylation[!is.na(total_counts_m_n), ]
      rm(df_phased_methylation, fls)
      
      # get allele-specific methylation
      df_tumour[, CN := nA+nB]
      df_phased_methylation <- df_tumour[order(chrom,start,end,width,POS,genotype), 
                                         .(genotype_1="ref", genotype_2="alt", num_genotypes=length(genotype),
                                           POS_genotype=POS_genotype[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                           total_depth_1 = total_depth[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                           total_depth_2 = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                           POS_genotype_BAF_t = total_depth[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_depth),
                                           total_counts_1 = total_counts[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                           total_counts_2 = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                           CpG_BAF_t = total_counts[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]/sum(total_counts),
                                           M_1 = M[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]], 
                                           UM_1 = UM[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                           total_counts_m_1 = total_counts_m[genotype=="ref"&POS_genotype%in%POS_genotype[genotype=="alt"]],
                                           M_2 = M[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]], 
                                           UM_2 = UM[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]],
                                           total_counts_m_2 = total_counts_m[genotype=="alt"&POS_genotype%in%POS_genotype[genotype=="ref"]]),
                                         keyby=c("CHR","chrom","start","end","width","POS","ref","alt", "nA", "nB", "CN")]  
      
      # label and order cols
      df_phased_methylation[, label := ifelse((POS==start & ref=="C")|(POS==end & ref=="G"), "CG-destroying","CG-forming")]
      df_phased_methylation <- 
        df_phased_methylation[, .SD, .SDcols=c(1:which(colnames(df_phased_methylation)=="CpG_BAF_t"), which(colnames(df_phased_methylation)=="label"),
                                               (which(colnames(df_phased_methylation)=="CpG_BAF_t")+1):which(colnames(df_phased_methylation)=="total_counts_m_2"))]
      
      # join tumour and normal info
      df_phased_methylation <- merge(df_phased_methylation, df_normal, 
                                     by = c("CHR", "chrom", "start", "end", "width","POS", "ref", "alt", 
                                            "POS_genotype", "genotype_1", "genotype_2", "num_genotypes"), all=FALSE)
      rm(df_normal)
      
      # get phased methylation rates
      df_phased_methylation[, m_1 := as.numeric(M_1)/as.numeric(total_counts_m_1)] 
      df_phased_methylation[, m_2 := as.numeric(M_2)/as.numeric(total_counts_m_2)]
      
      # extract purity
      rho <- metadata$purity[jobid]
      rho <- as.numeric(unname(unlist(rho)))
      
      # # load phased haplotypes
      # fls = paste0(parent_dir,PT_ID,"/Copy_number/",PT_ID,"_haplotypes")
      # fls = paste0(fls,".RData")
      # load(fls)
      
      # # set for overlap 
      # dt_phasing <- data.table(dt_phasing)
      # dt_phasing[, c("POS_genotype") := list(as.character(POS))]
      # dt_phasing <- dt_phasing[, .SD, .SDcols=c("chrom", "POS_genotype", paste0("phased_",sample_id), "phasing", "re_seg")]
      # colnames(dt_phasing)[colnames(dt_phasing)==paste0("phased_",sample_id)] <- "phased_seg_BAF"
      
      # # overlap phased methylation rates with haplotypes
      # df_phased_methylation[, POS_genotype := tstrsplit(POS_genotype, ",", keep=1)]
      # tmp <- merge(df_phased_methylation, dt_phasing, by = c("chrom", "POS_genotype"), all=FALSE)
      # rm(dt_phasing)



      ### this is added by lesky: if there is no haplotype object, just add a column "phasing" to be "unknown", test the remaining code
      tmp <- df_phased_methylation
      tmp$phasing <- "unknown"
      ### end of added lines
      
      # get purified tumour phased methylation rates
      df_phased_methylation <- get_phased_pure_tumour_methylome(tmp, p=rho, sex = sex)
      rm(tmp)
      
      # Checkpoint
      print("Purified tumour methylation rate HDI added for each CpG locus per parental alleles")
      
      # bound values between 0 and 1
      df_phased_methylation[, m_t_1_raw := m_t_1] 
      df_phased_methylation[, m_t_2_raw := m_t_2]
      df_phased_methylation[, m_t_1 := ifelse(m_t_1_raw<0, 0, 
                                              ifelse(m_t_1_raw>1, 1, m_t_1_raw))] 
      df_phased_methylation[, m_t_2 := ifelse(m_t_2_raw<0, 0, 
                                              ifelse(m_t_2_raw>1, 1, m_t_2_raw))]
      
      # get potentially imprinted loci
      df_phased_methylation[,m_diff_t12 := abs(m_t_1-m_t_2)]
      df_phased_methylation[,m_diff_n12 := abs(m_n_1-m_n_2)]
      df_phased_methylation[,m_diff_b12 := abs(m_1-m_2)]
      
      # get methylation difference
      df_phased_methylation[,m_diff_tn1 := abs(m_n_1-m_t_1)]
      df_phased_methylation[,m_diff_tn2 := abs(m_n_2-m_t_2)]
      df_phased_methylation[,m_diff_bn1 := abs(m_n_1-m_1)]
      df_phased_methylation[,m_diff_bn2 := abs(m_n_2-m_2)]
      
      # Compute the tumour methylation rate HDI 
      n <- nrow(df_phased_methylation)
      df_phased_methylation[, as.character(c("m_t_1_low", "m_t_1_high")) := 
                              as.list(HDIofMCMC(M_b=M_1,UM_b=UM_1,M_n=M_n_1,UM_n=UM_n_1,p=rho,CN=n1,CN_n=n1_n,0.99)), by=.(1:n)]
      
      # cap allele 1 HDI values
      df_phased_methylation[, m_t_1_low := ifelse(m_t_1_low < 0, 0, ifelse(m_t_1_low > 1 ,1, m_t_1_low))]
      df_phased_methylation[, m_t_1_high := ifelse(m_t_1_high < 0, 0, ifelse(m_t_1_high > 1 ,1, m_t_1_high))]
      
      # Compute the tumour methylation rate HDI 
      ## Allele 2
      df_phased_methylation[, as.character(c("m_t_2_low", "m_t_2_high")) := 
                              as.list(HDIofMCMC(M_b=M_2,UM_b=UM_2,M_n=M_n_2,UM_n=UM_n_2,p=rho,CN=n2,CN_n=n2_n,0.99)), by=.(1:n)]
      rm(n)
      
      # cap allele 2 HDI values
      df_phased_methylation[, m_t_2_low := ifelse(m_t_2_low < 0, 0, ifelse(m_t_2_low > 1 ,1, m_t_2_low))]
      df_phased_methylation[, m_t_2_high := ifelse(m_t_2_high < 0, 0, ifelse(m_t_2_high > 1 ,1, m_t_2_high))]
      
      # Add allele 1/2 DMPs
      df_phased_methylation[, m_t_2_high := ifelse(m_t_2_high < 0, 0, ifelse(m_t_2_high > 1 ,1, m_t_2_high))]
      
      # check overlaps
      df_phased_methylation[, is_DM_t12 :=  ((m_t_1_low>m_t_2_high)|(m_t_2_low>m_t_1_high)) & m_diff_t12 > 0.5]
      #df_phased_methylation[!is.na(is_DM_t12)&ov_t12<0.21&m_diff_t12>0.5, table(is_DM_t12==TRUE,round2(ov_t12,digits=2))] 
      
      # Checkpoint
      print("DM between alleles obtained using HDI values")
      
      # join phased and bulk tumour info
      df_phased_methylation <- merge(df_phased_methylation, dt2, by = c("chrom", "start", "end"), all=FALSE)
      
      # order colnames
      df_phased_methylation <- 
        df_phased_methylation[, .SD, .SDcols=c(1:which(colnames(df_phased_methylation)=="label"),
                                               which(colnames(df_phased_methylation)=="total_counts_n_1"):which(colnames(df_phased_methylation)=="CpG_BAF_n"),
                                               which(colnames(df_phased_methylation)=="n1"):which(colnames(df_phased_methylation)=="n2_n"),
                                               which(colnames(df_phased_methylation)=="M_1"):which(colnames(df_phased_methylation)=="total_counts_m_2"),
                                               which(colnames(df_phased_methylation)=="m_1"):which(colnames(df_phased_methylation)=="m_2"),
                                               which(colnames(df_phased_methylation)=="m_b"):which(colnames(df_phased_methylation)=="m_b_high"),
                                               which(colnames(df_phased_methylation)=="M_n_1"):which(colnames(df_phased_methylation)=="m_n_high"), 
                                               which(colnames(df_phased_methylation)=="m_t_1"):which(colnames(df_phased_methylation)=="is_DM_t12"), 
                                               which(colnames(df_phased_methylation)=="m_t_corr"):which(colnames(df_phased_methylation)=="DMP_t"))]
      # df_phased_methylation[!is.na(m_t_2)&!is.na(m_t_1),table(m_diff_t12 > 0.6, !is.na(DMP_t))]
      # df_phased_methylation[!is.na(m_t_2)&!is.na(m_t_1),table((m_t_1_low>m_t_2_high)|(m_t_2_low>m_t_1_high), !is.na(DMP_t))]
      # df_phased_methylation[!is.na(m_t_2)&!is.na(m_t_1),table(!is.na(is_DM_t12)&is_DM_t12==TRUE, !is.na(DMP_t))]
      
      
      #save
      save(df_phased_methylation, file=paste0("/rsrch4/home/genetics/slai4/CAMDAC_paper_revision/ASM_imprinting/ASM_RRBS/", PT_ID,
                                                "/T/allele_specific_methylation_rates_T_v1.RData"))
    }
  }
})
}    
