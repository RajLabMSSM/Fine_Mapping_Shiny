
library(shiny)
library(dplyr)
# devtools::install_github('rstudio/DT')
# library(echolocatoR)
library(DT)
library(data.table)
library(ggplot2)
library(plotly) 


printer <- function(..., v=T){if(v){print(paste(...))}}

construct_SNPs_labels <- function(subset_DT, lead=T, method=T, consensus=T, verbose=F, remove_duplicates=T){ 
  printer("+ PLOT:: Constructing SNP labels...", v=verbose)
  labelSNPs <- data.table::data.table() 
  subset_DT <- data.table::as.data.table(subset_DT)
  CS_cols <- colnames(subset_DT)[endsWith(colnames(subset_DT),".Credible_Set")]
  methods <- gsub("\\.Credible_Set","",CS_cols)
  
  ## BEFORE fine-mapping  
  if(lead){
    before <- subset( subset_DT %>% arrange(P), leadSNP == T)
    before$type <- "Lead SNP"
    before$color <- "red"
    before$shape <- 18
    before$size <- 2
    labelSNPs <- rbind(labelSNPs, before, fill=T)
  }
  if(method){
    # AFTER fine-mapping
    after = subset(subset_DT, Support>0) 
    if(dim(after)[1]>0){
      after$type <- "Credible Set"  
      after$color<- "green3"
      after$shape <- 1
      after$size=4 
      after$Credible_Sets <- lapply(1:nrow(after), function(i){
        paste(methods[data.frame(after)[i,CS_cols]>0], collapse=", ")
      }) %>% unlist() 
      labelSNPs <- rbind(labelSNPs, after, fill=T) 
    } 
  } 
  if(consensus & "Consensus_SNP" %in% colnames(subset_DT)){
    # Conensus across all fine-mapping tools
    cons_SNPs <- subset(subset_DT, Consensus_SNP==T)
    if(dim(cons_SNPs)[1]>0){
      cons_SNPs$type <- "Consensus SNP"
      cons_SNPs$color <- "darkgoldenrod1"
      cons_SNPs$shape <- 1
      cons_SNPs$size=6
      cons_SNPs$Credible_Sets <- lapply(1:nrow(cons_SNPs), function(i){
        paste(methods[data.frame(cons_SNPs)[i,CS_cols]>0], collapse=", ")
      }) %>% unlist() 
      labelSNPs <- rbind(labelSNPs, cons_SNPs, fill=T)
    } 
  } 
  # Convert to GRanges object
  # labelGR <- transformDfToGr(data=labelSNPs, seqnames = "CHR", start = "POS", end = "POS")
  # names(labelGR) <- labelGR$SNP
  # plotGrandLinear(gr.snp, aes(y = P, x=POS), highlight.gr = labelGR)
  
  # If there's duplicates only show the last one
  if(remove_duplicates){
    labelSNPs$rowID <- 1:nrow(labelSNPs)
    labelSNPs <- labelSNPs %>% 
      dplyr::group_by(SNP) %>% 
      dplyr::arrange(rowID) %>% 
      dplyr::slice(n())  
  } else {
    labelSNPs$rowID <- 1:nrow(labelSNPs)
    labelSNPs <- labelSNPs %>% 
      dplyr::group_by(SNP, type) %>% 
      dplyr::arrange(rowID) %>% 
      dplyr::slice(n())  
  }
  
  labelSNPs$type <- factor(x = labelSNPs$type, 
                           levels = c("Lead SNP","Credible Set","Consensus SNP"), 
                           ordered = T)
  return(as.data.frame(labelSNPs))
}


find_consensus_SNPs <- function(finemap_DT,
                                verbose=T,
                                credset_thresh=.95,
                                consensus_thresh=2,
                                sort_by_support=T,
                                exclude_methods=NULL){
  printer("+ Identifying Consensus SNPs...",v=verbose)
  exclude_methods <- append(exclude_methods,"mean")
  # Find SNPs that are in the credible set for all fine-mapping tools
  CS_cols <- colnames(finemap_DT)[endsWith(colnames(finemap_DT),".Credible_Set")]
  CS_cols <- CS_cols[!(CS_cols %in% paste0(exclude_methods,".Credible_Set"))]
  if(consensus_thresh=="all"){consensus_thresh<-length(CS_cols)}
  printer("++ support_thresh =",consensus_thresh)
  # Get the number of tools supporting each SNP
  ## Make sure each CS is set to 1
  support_sub <- subset(finemap_DT, select = CS_cols) %>% data.frame()
  support_sub[sapply(support_sub, function(e){e>1})] <- 1
  finemap_DT$Support <- rowSums(support_sub, na.rm = T)
  finemap_DT$Consensus_SNP <- finemap_DT$Support >= consensus_thresh
  # Sort
  if(sort_by_support){
    finemap_DT <- finemap_DT %>% arrange(desc(Consensus_SNP), desc(Support))
  }
  
  # Calculate mean PP
  printer("+ Calculating mean Posterior Probability (mean.PP)...")
  PP.cols <- grep(".PP",colnames(finemap_DT), value = T)
  PP.cols <- PP.cols[!(PP.cols %in% paste0(exclude_methods,".PP"))]
  PP.sub <- subset(finemap_DT, select=c("SNP",PP.cols)) %>% data.frame()# %>% unique()
  PP.sub[is.na(PP.sub)] <- 0
  if(NCOL(PP.sub[,-1]) > 1){
    finemap_DT$mean.PP <- rowMeans(PP.sub[,-1], na.rm = T)
  } else{
    finemap_DT$mean.PP <- PP.sub[,-1]
  }
  finemap_DT$mean.Credible_Set <- ifelse(finemap_DT$mean.PP>=credset_thresh,1,0)
  
  # PP.sub %>% arrange(desc(mean.PP)) %>% head()
  printer("++",length(CS_cols),"fine-mapping methods used.")
  printer("++",dim(subset(finemap_DT,Support>0))[1],"Credible Set SNPs identified.")
  printer("++",dim(subset(finemap_DT,Consensus_SNP==T))[1],"Consensus SNPs identified.")
  return(finemap_DT)
}



assign_lead_SNP <- function(new_DT, verbose=T){
  if(sum(new_DT$leadSNP)==0){
    printer("+ leadSNP missing. Assigning new one by min p-value.", v=verbose)
    top.snp <- head(arrange(new_DT, P, desc(Effect)))[1,]$SNP
    new_DT$leadSNP <- ifelse(new_DT$SNP==top.snp,T,F)
  }
  return(new_DT)
}

make_locus_dict <- function(root="~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019",
                            pattern="*_ggbio*.png|multi_finemap_plot.png",
                            slice_n=1){ 
  locus_plots_df <- data.frame(path=list.files(path = root, 
                                               pattern = pattern,
                                               full.names = T, recursive = T), stringsAsFactors = F) %>% 
    dplyr::mutate(locus=basename(dirname(dirname(path)))) %>%
    dplyr::arrange(path) %>%
    dplyr::group_by(locus) %>% 
    dplyr::slice(slice_n) 
  locus_plots <- setNames(locus_plots_df$path, locus_plots_df$locus)
  return(locus_plots)
}


# Only need to run this once beforehand
prepare_files <- function(root="~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"){
  
  # ------- Collect plot paths -------
  locus_plots <- make_locus_dict(root=root,
                                 pattern="_ggbio.png")#multi_finemap_plot.png  
  dir.create("www/plots", showWarnings = F, recursive = T)
  locus_plots <- lapply(names(locus_plots), function(locus){
    new_path <- file.path("www","plots",paste0(locus,"_plot.png"))
    file.copy(from = locus_plots[[locus]], 
              to = new_path, overwrite = T)
    return(new_path)
  }) %>% unlist() %>% `names<-`(names(locus_plots))
  saveRDS(locus_plots, "www/locus_plots.RDS")
  
  
  # ------- Collect table paths -------
  locus_tables <- make_locus_dict(root=root,
                                  pattern="Multi-finemap_results.txt")   
  leadSNPs <<- c()
  dir.create("www/data", showWarnings = F, recursive = T)
  locus_tables <- lapply(names(locus_tables), function(locus){
    message(locus)
    new_path <- file.path("www","data",paste0(locus,"_data.csv"))
    dat <- data.table::fread(locus_tables[[locus]])
    
    # Do preprocessing beforehand 
    if(!"Locus" %in% colnames(dat)){
      dat <- cbind(Locus=locus, dat)
    }
    dat$Mb <- dat$POS/1000000  
    dat <- assign_lead_SNP(dat)
    dat <- find_consensus_SNPs(dat, verbose = F) 
    dat$proportion_cases <- round(dat$proportion_cases, 5)
    dat$mean.PP <- round(dat$mean.PP, 5)
    
    write.csv(dat, new_path, row.names = F)
    # data.table::fwrite(data.table::data.table(dat), new_path,) 
    
    # Store lead SNP info
    leadSNPs <<- append(leadSNPs, subset(dat, leadSNP)$SNP[1])
    return(new_path)
  }) %>% unlist() %>% `names<-`(names(locus_tables))
  
  names(leadSNPs) <-  names(locus_tables)
  saveRDS(locus_tables, "www/locus_tables.RDS")
  
  
  
  # -------Collect LD paths -------
  locus_LD <- make_locus_dict(root=root,
                              pattern="UKB_LD.RDS" )
  locus_LD <- locus_LD[names(locus_LD)!="plink"]
  dir.create("www/LD", showWarnings = F, recursive = T)
  locus_LD <- parallel::mclapply(names(locus_LD), function(locus){
    message(locus)
    new_path <- file.path("www","LD",paste0(locus,"_LD.csv"))
    # Get lead SNP that's ALSO in LD_matrix
    # dat <- data.table::fread(locus_tables[[locus]]) 
    # dat <- assign_lead_SNP(subset(dat, SNP %in% colnames(LD_matrix)))
    # write.csv(dat, locus_tables[[locus]]) # rwrite with new lead SNP
    # lead_snp <- subset(dat, leadSNP)$SNP[1]
    if(endsWith(locus_LD[[locus]], suffix = ".RDS")){
      LD_matrix <- readRDS(locus_LD[[locus]])
    }
    if(endsWith(locus_LD[[locus]], suffix = ".RData")){
      load(locus_LD[[locus]])
    }
    
    LD_df <-  tryCatch(expr = { 
      LD_df <- data.frame(LD_matrix[,leadSNPs[[locus]] ])
      LD_df <- cbind(SNP=colnames(LD_matrix), LD_df)
      colnames(LD_df)[2] <- leadSNPs[[locus]]
      LD_df
    }, 
    error = function(e){
      LD_df <- data.frame(SNP=colnames(LD_matrix),
                          r=rep(NA,length( colnames(LD_matrix))))
      colnames(LD_df)[2] <- leadSNPs[[locus]]  
      LD_df
    }, finally = {
      write.csv(LD_df, new_path, row.names = F) 
      LD_df 
    }) 
    return(new_path)
  }, mc.cores = 4) %>% unlist() %>% `names<-`(names(locus_LD))
  saveRDS(locus_LD, "www/locus_LD.RDS") 
}


locus_plots <- readRDS("www/locus_plots.RDS") 
locus_tables <- readRDS("www/locus_tables.RDS")
locus_LD <- readRDS("www/locus_LD.RDS")


# Dropdown bars
studies <- "Nalls23andMe_2019" #list.dirs("./Data/GWAS/", full.names = F, recursive = F)
loci <-  sort(setNames(names(locus_tables), names(locus_tables) ) )

#### --------- UI ---------- ####
ui <- fluidPage(
  
  sidebarLayout(position = "left", 
                sidebarPanel(style = "position:fixed; width:25%; height:95vh; overflow-y:auto",
                             h2("Fine-mapping Results Browser"),
                             hr(), 
                             selectInput(inputId = "study", 
                                         label = "Study : ", 
                                         choices = studies, 
                                         selected = "Nalls23andMe_2019"),
                             selectInput(inputId = "locus", 
                                         label = "Locus : ",  
                                         choices = loci, 
                                         selected = "LRRK2", #loci[[1]]
                             ),
                             shiny::checkboxInput(inputId = "separate_finemap_methods", 
                                                  label = "Separate fine-mapping methods", 
                                                  value = FALSE),
                             hr(),
                             h4("Interactive plot directions"),
                             shiny::p("Layered results from fine-mapping each locus."), 
                             shiny::p("Hover over each point (SNP) to show more details."),
                             shiny::p("Click and drag to zoom in. Double click to zoom out."),
                             h5("Key:"),
                             span(strong("r2", style="display: inline;"),p("Pairwise LD correlation with lead GWAS SNP", style="display: inline;")),br(),
                             span(p("◆",style="color:red; font-size:25px; display: inline"), p("Lead GWAS SNP", style="display: inline;")),br(),
                             span(p("○",style="color:green; font-size:25px; display: inline"), p("Union Credible Set SNP", style="display: inline;")),br(),
                             span(p("○",style="color:goldenrod; font-size:25px; display: inline"), p("Consensus SNP", style="display: inline;")),br(),
                             # tags$ul(
                             #   tags$li(span(p("◆",style="color:red, size:10"), p("Lead GWAS SNP"))),
                             #   tags$li(span(p("○",style="color:green"), p("Union Credible Set SNP")),
                             #   tags$li(p("○",style="color:goldenrod"), p("Consensus SNP")),
                             # ),
                             h4("Fine-mapping results"),
                             shiny::p("Raw summary statistics and fine-mapping results for each locus."),
                             
                             hr(),
                             h4("Authors"),
                             p("Raj Lab"),
                             p("Dept. of Neuroscience, Dept. of Genetics & Genomic Sciences"),
                             p("Icahn School of Medicine at Mount Sinai"),
                             p("New York, NY")
                ),
                mainPanel( 
                  fluidRow(
                    column(width = 12,
                           h4("Interactive plot"),   
                           h5(textOutput("locus_name"), align="center"),
                           h6(textOutput("n_snps"), align="center"),
                           # plotOutput("plot")
                           plotly::plotlyOutput("plotly", 
                                                height = "500px")
                    ) 
                  ), 
                  fluidRow( 
                    column(width = 12, 
                           h4("Static plot"), 
                           plotOutput("plot", inline = T) 
                    )
                  ),
                  br(),
                  fluidRow( 
                    column(width = 12, 
                           h4("Fine-mapping results"), 
                           DT::dataTableOutput("results")
                           # tableOutput("results")
                    ) 
                  ),
                  br(),
                ),
                
  ),
)




#### ---------SERVER ####
server <- function(input, output, session) { 
  
  import_data <- function(input){
    # finemap_DT <- data.table::fread("./www/data/WNT3_data.csv")
    finemap_DT <- data.table::fread(locus_tables[[input$locus]] ) 
    # Import and merged LD info
    LD_df <- data.table::fread(locus_LD[[input$locus]]) 
    LD_df$r2 <- LD_df[,2]^2 
    finemap_DT <- data.table::merge.data.table(x = finemap_DT,
                                               y = LD_df,
                                               by = "SNP")
    return(finemap_DT)
  }
  
  # Plot  
  output$plot <- renderImage({
    # print(locus_plots[[input$locus]]) 
    plot_info <- tryCatch(expr = {list(src = locus_plots[[input$locus]], width="100%") },
                          error=function(e){
                            list(src = "./shiny_input/icons/under_construction.png", width="40px")
                          }) 
    return(plot_info)
  }, deleteFile = FALSE)
  
  output$locus_name <- renderText({input$locus})
  
  output$plotly <- plotly::renderPlotly({
    withProgress(message = 'Loading', value = 0, {
      finemap_DT <- import_data(input)
      finemap_methods <- gsub("\\.PP", "",grep(pattern = "\\.PP", colnames(finemap_DT), value = T ))
      output$n_snps <- renderText({ paste(length(unique(finemap_DT$SNP)),"SNPs") })
      
      
      n_steps=1+1+length(finemap_methods)+1+1
      incProgress(1/n_steps, detail = paste("Importing data...")) 
      
      
      snp_plot <- function(finemap_DT, 
                           y_var="-log10(P)",
                           ylab_prefix=NULL,
                           locus=NULL, 
                           interactive=T,
                           viridis_color=F,
                           ylimits=NULL){
        snp.labs <- construct_SNPs_labels(finemap_DT,
                                          remove_duplicates = F)
        
        gp <- ggplot(data=finemap_DT, aes_string(x="Mb", y=y_var, 
                                                 color="r2",
                                                 label_CHR="CHR",
                                                 label_SNP="SNP",
                                                 label_Effect="Effect",
                                                 label_P="P", 
                                                 label_StdErr="StdErr",
                                                 label_A1="A1", 
                                                 label_A2="A2") ) + 
          geom_point(alpha=.5, show.legend = T, size=2) +
          geom_point(data = snp.labs,
                     aes_string(x="Mb", y=y_var,
                                label_type="type", 
                                label_CS="Credible_Sets"), 
                     shape=snp.labs$shape, 
                     color=snp.labs$color, 
                     size=snp.labs$size, 
                     alpha=.7, show.legend = F) + 
          labs(title =locus,# ifelse(is.null(locus), NULL, paste("Locus :", locus)), 
               y=paste(ylab_prefix,y_var)) + 
          scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0,.2)))  +
          scale_color_gradient(low = "blue", high = "red", na.value = "gray", 
                               limits=c(0,1), breaks=c(0,.25,.5,.75,1)) +
          theme_bw() +
          theme(axis.title.y = element_text(size=8), 
                plot.title = element_text(hjust = .5)) 
        
        if(!is.null(ylimits)){
          gp <- gp + ylim(ylimits)
        }
        if(viridis_color){
          gp <- gp + scale_color_viridis_c()
        }  
        if(interactive){
          return(plotly::ggplotly(gp))
        } else {return(gp)}
      }
      
      plt_list <- c()
      incProgress(1/n_steps, detail = paste("Creating GWAS plot..."))
      plt_list[["gwas"]] <- snp_plot(finemap_DT, 
                                     y_var="-log10(P)", 
                                     ylab_prefix = "GWAS")  
      if(input$separate_finemap_methods){
        for(m in finemap_methods){
          incProgress(1/n_steps, detail = paste("Creating",m,"plot..."))
          plt_list[[m]] <- snp_plot(finemap_DT, 
                                    y_var=paste0(m,".PP"), 
                                    viridis_color = F,
                                    ylimits = c(0,1.1))
        }  
      }  
      
      # output$plotly_height <- renderUI({
      #   if(input$separate_finemap_methods){"600px"}else{"400px"}
      # })
      
      pltly <- plotly::subplot(plt_list, 
                               nrows = length(plt_list), 
                               shareY = F, shareX = T,
                               titleY = T)
      incProgress(1/n_steps, detail = paste("Rendering merged plots...")) 
      return(pltly) 
      # incProgress(1/n_steps, detail = paste("Complete!")) 
    }) 
  }) # withProgress
  
  
  # TABLE 
  getOpts <- function(file_name="finemapping_results", 
                      rowGroup=NULL){
    opts <- list(scrollY = 500, 
                 sScrollX="100%", 
                 scrollX=T,
                 bScrollCollapse=T,
                 pageLength=50,
                 paging=F,
                 dom = 'frtipB', 
                 buttons = list( list(extend = 'csv', filename=file_name),
                                 list(extend = 'excel', filename=file_name),
                                 list(extend = 'pdf', filename=file_name),
                                 list(extend = 'print', filename=file_name),
                                 'copy'),  
                 # 1st column is 1 ONLY if rownames=F
                 fixedColumns = list(leftColumns = 2),
                 rowGroup = rowGroup)
    return(opts)
  }
  
  # Results table
  output$results <- DT::renderDT({  
    finemap_DT <- import_data(input)
    # Columns are zero-indexed
    consensus_index <- grep("Consensus_SNP", colnames(finemap_DT))-1
    
    DT::datatable(data = finemap_DT, 
                  extensions = c("Buttons","Scroller","FixedColumns","FixedHeader","RowGroup"),
                  class='cell-border stripe compact table-hover',
                  options = getOpts(file_name = paste0(input$locus,"results",sep="_"),
                                    # Can't have rowGroup and FixedColumns at the same time
                                    rowGroup = NULL#list(dataSrc = consensus_index)
                  ),
                  filter='top', 
                  rownames=F, 
                  selection='single')
  })
}


shinyApp(ui, server)


