
# ---- Functions for preprocessing echolocatoR output ---#
## Handles much of the ambiguity in file naming schemes,
## since echolocatoR has changed over time.
## Also handles results from multiple LD panels.

library(dplyr)
# library(data.table)
# library(echolocatoR)


prepare_files <- function(root="/sc/arion/projects/pd-omics/brian/Fine_Mapping"){
  ####  Collect plot paths #### 
  locus_plots_df <- prepare_plots(root=root, 
                                  pattern="*_ggbio*.png|multi_finemap_plot.png|^multiview\\.",
                                  force_new_plot = T) 
  #### Collect table paths #### 
  locus_tables_df <- prepare_tables(root=root,
                                    pattern="*\\.Multi-finemap.tsv.gz",
                                    force_new_table = T)  
  #### Collect LD paths ####
  locus_LD_df <- prepare_LD(root=root,
                            locus_tables_df=locus_tables_df,
                            pattern="*UKB_LD.RDS|*1KGphase3_LD.RDS|*1KGphase1_LD.RDS",
                            force_new_ld=F)
  #### Gather processed files ####
  all_paths <- gather_processed_paths(processed_dir="www/data") 
  saveRDS(all_paths, "www/all_paths.RDS")
  return(all_paths)
}




infer_LD_panel <- function(path){
  if(grepl("1KGphase3", path)) return("1KGphase3")
  if(grepl("1KGphase1", path)) return("1KGphase1")
  if(grepl("UKB", path)) return("UKB")
  return(NA)
}


make_locus_df <- function(root="/sc/arion/projects/pd-omics/brian/Fine_Mapping",
                          pattern,
                          slice_n=NULL){ 
  locus_df <- data.frame(path=list.files(path = root, 
                                               pattern = pattern,
                                               full.names = T, recursive = T), stringsAsFactors = F) %>% 
    dplyr::mutate(subfolder=basename(dirname(path)) %in% c("Multi-finemap","LD","plink") ) %>% 
    # Some plots are in their own subfolder, others are not
    dplyr::mutate(locus=ifelse(subfolder, basename(dirname(dirname(path))), basename(dirname(path))),
                  dataset=ifelse(subfolder, basename(dirname(dirname(dirname(path)))), basename(dirname(dirname(path))) ),
                  dataset_type=ifelse(subfolder, basename(dirname(dirname(dirname(dirname(path))))), basename(dirname(dirname(dirname(path)))) ) 
    )
   
  locus_df$LD_ref <- lapply(basename(locus_df$path), function(x){infer_LD_panel(x)}) %>% unlist() 
  locus_df <- locus_df %>%  tidyr::separate(path, sep = "[.]", into = c(NA,NA,NA,"zoom"), remove = F)
  
  locus_df <- locus_df %>% 
    dplyr::mutate(locus_dir=file.path("www/data",dataset_type,dataset,locus)) %>%
    dplyr::mutate(plot_path=file.path(locus_dir,"plots",paste(locus,LD_ref,"locus_plot",zoom,"png",sep=".")),
                  data_path=file.path(locus_dir,"multi_finemap",paste(locus,LD_ref,"multi_finemap","csv.gz",sep=".")),
                  ld_path=file.path(locus_dir,"LD",paste(locus,LD_ref,"LD","csv.gz",sep="."))
    ) 
  # Arbitrarily use one plot per Locus
  if(!is.null(slice_n)){
    locus_df <- locus_df %>%
      dplyr::group_by(locus, dataset_type, LD_ref, zoom) %>% 
      dplyr::slice(slice_n) 
  }   
  return(data.table::data.table(locus_df))
}


prepare_plots <- function(root,
                          pattern="*_ggbio*.png|multi_finemap_plot.png|^multiview\\.",
                          force_new_plot = T){
  locus_plots_df <- make_locus_df(root=root, 
                                  pattern=pattern)  
  new_plots <- lapply(1:nrow(locus_plots_df), function(i){
    ROW <- locus_plots_df[i,] 
    print(paste(basename(ROW$path),"==>",basename(ROW$plot_path)))
    dir.create(dirname(ROW$plot_path), showWarnings = F, recursive = T)
    if((!file.exists(ROW$plot_path)) | force_new_plot ){
      file.copy(from = ROW$path, 
                to = ROW$plot_path, 
                overwrite = T)
    }   
    return(ROW$plot_path)
  }) %>% unlist() 
  return(locus_plots_df)
} 



prepare_tables <- function(root,
                           pattern="*\\.Multi-finemap.tsv.gz",
                           force_new_table = T){
  locus_tables_df <- make_locus_df(root=root,
                                   pattern=pattern)   
  # !!!! requires echolocatoR  !!!! 
  locus_tables_df$leadSNP <- parallel::mclapply(1:nrow(locus_tables_df), function(i){
    ROW <- locus_tables_df[i,]
    print(paste(basename(ROW$path),"==>",basename(ROW$data_path)))
    if((!file.exists(ROW$path)) | force_new_table){ 
      dat <- data.table::fread(ROW$path, nThread = 1)
      # Do preprocessing beforehand 
      if(!"Locus" %in% colnames(dat)){
        dat <- cbind(Locus=ROW$locus, dat)
      }
      dat <- echolocatoR::update_cols(finemap_dat = dat)
      dat$Mb <- dat$POS/1000000  
      dat <- echolocatoR::assign_lead_SNP(dat)
      dat <- echolocatoR::find_consensus_SNPs(dat, verbose = F) 
      if("proportion_cases" %in% colnames(dat)){
        dat$proportion_cases <- round(dat$proportion_cases, 5)
      }
      dat$mean.PP <- round(dat$mean.PP, 5)
      # Write standardized file
      dir.create(dirname(ROW$data_path), showWarnings = F, recursive = T)
      data.table::fwrite(dat, ROW$data_path, row.names = F, nThread = 1)  
      leadSNP <- subset(dat, leadSNP)$SNP[1]
      return(leadSNP)
    } else{
      dat <- data.table::fread(ROW$data_path, nThread = 1)
      leadSNP <- subset(dat, leadSNP)$SNP[1]
      return(leadSNP)
    }
  }, mc.cores = 1 #parallel::detectCores()
  ) %>% unlist() 
  return(locus_tables_df)
}


prepare_LD <- function(root,
                       locus_tables_df,
                       pattern="*UKB_LD.RDS|*1KGphase3_LD.RDS|*1KGphase1_LD.RDS",
                       force_new_ld=F){
  locus_LD_df <- make_locus_df(root=root,
                               pattern=pattern) 
  ld_paths <- parallel::mclapply(1:nrow(locus_LD_df), function(i,
                                                               .locus_tables_df=locus_tables_df){
    ROW <- locus_LD_df[i,]
    locus <- ROW$locus
    print(paste(basename(ROW$path),"==>",basename(ROW$ld_path)))
    if((!file.exists(ROW$ld_path)) | force_new_ld){
      LD_df <- tryCatch(expr = {
        # Get lead SNP that's ALSO in LD_matrix 
        if(endsWith(ROW$path, suffix = ".RDS")){
          LD_matrix <- readRDS(ROW$path)
        }
        if(endsWith(ROW$path, suffix = ".RData")){
          load(ROW$path)
        } 
        # Get the lead SNP for this dataset's locus (ID'ed by locus_dir)
        tryCatch({
          lead_snp <- subset(.locus_tables_df, locus_dir==ROW$locus_dir)$leadSNP 
        })
        
        LD_df <-  tryCatch(expr = {  
          LD_df <- data.frame(LD_matrix[, lead_snp])
          LD_df <- cbind(SNP=colnames(LD_matrix), LD_df)
          colnames(LD_df)[2] <- lead_snp
          return(LD_df)
        }, 
        error = function(e){
          LD_df <- tryCatch(expr = {
            LD_df <- data.frame(SNP=colnames(LD_matrix),
                                r=rep(NA,length( colnames(LD_matrix))))
            colnames(LD_df)[2] <-lead_snp
            return(LD_df)
          }, 
          error = function(e){message("xxx LD failed @ processing xxx"); NULL})
        }, finally = {
          if(exists("LD_df")){
            dir.create(dirname(ROW$ld_path), showWarnings = F, recursive = T) 
            data.table::fwrite(LD_df, ROW$ld_path, row.names = F, nThread = 1) 
            return(LD_df) 
          } 
        })  
      },
      error =  function(e){message("xxx LD failed @ read in xxx"); NULL }) 
    } 
    return(ROW$ld_path)
  }, mc.cores = 1#parallel::detectCores() 
  ) %>% unlist() 
  return(locus_LD_df)
} 


gather_processed_paths <- function(processed_dir="www/data"){
  all_paths <- data.frame(file_path=list.files(path = processed_dir, full.names = T, recursive = T), stringsAsFactors = F) %>%
    dplyr::mutate(study_type=basename(dirname(dirname(dirname(dirname(file_path))))),
                  study=basename(dirname(dirname(dirname(file_path)))),
                  locus_dir=dirname(dirname(file_path)),
                  locus=basename(dirname(dirname(file_path))),
                  file_type=basename(dirname(file_path))) %>%
    # dplyr::mutate(zoom=ifelse(file_type=="plots", strsplit(basename(file_path),"\\.")[[1]][4], NA)
    #               ) %>%
    data.table::data.table()
  all_paths$LD_ref <- lapply(basename(all_paths$file_path), function(x){infer_LD_panel(x)}) %>% unlist()
  all_paths <- all_paths %>% tidyr::separate(file_path, sep = "[.]", into = c(NA,NA,NA,"zoom"), remove = F)
  all_paths[all_paths$file_type!="plots","zoom"] <- NA
  return(all_paths)
}

