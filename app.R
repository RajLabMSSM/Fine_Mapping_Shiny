
library(shiny)
library(shinyBS)
library(dplyr)
# devtools::install_github('rstudio/DT')
# library(echolocatoR)
library(DT)
library(data.table)
library(ggplot2)
library(plotly) 


get_study_metadata <- function(metadata_path="./www/metadata/GWAS-QTL_data_dictionary.xlsx"){
  gwas <- suppressMessages(readxl::read_excel(metadata_path, sheet = "GWAS"))
  qtl <- suppressMessages(readxl::read_excel(metadata_path, sheet = "QTL"))
  meta <- data.table::rbindlist(list(cbind(dataset_type="GWAS", gwas), 
                                     cbind(dataset_type="QTL",qtl)), fill=T) %>%
    dplyr::select(dataset_type, dataset, phenotype, prop_cases, build, reference)
  return(meta) 
}

label_studies <- function(all_paths){
  meta <- get_study_metadata()
  meta_sub <- subset(meta, dataset %in% gsub("\\.n_causal1","",unique(all_paths$study)))
  return(meta_sub)
}

make_locus_df <- function(root="/Volumes/Scizor/Fine_Mapping/Data",
                          pattern="*_ggbio*.png|multi_finemap_plot.png|^multiview\\.",
                          slice_n=1){ 
  locus_plots_df <- data.frame(path=list.files(path = root, 
                                               pattern = pattern,
                                               full.names = T, recursive = T), stringsAsFactors = F) %>% 
    dplyr::mutate(subfolder=basename(dirname(path)) %in% c("Multi-finemap","LD","plink") ) %>% 
    # Some plots are in their own subfolder, others are not
    dplyr::mutate(locus=ifelse(subfolder, basename(dirname(dirname(path))), basename(dirname(path))),
                  dataset=ifelse(subfolder, basename(dirname(dirname(dirname(path)))), basename(dirname(dirname(path))) ),
                  dataset_type=ifelse(subfolder, basename(dirname(dirname(dirname(dirname(path))))), basename(dirname(dirname(dirname(path)))) ) 
                  )
      
  locus_plots_df$LD_ref <- ifelse(startsWith("multiview\\.",basename(locus_plots_df$path)), strsplit( basename(locus_plots_df$path), "\\.")[[1]][3], NA)
  locus_plots_df$zoom <- ifelse(startsWith("multiview\\.",basename(locus_plots_df$path)), strsplit( basename(locus_plots_df$path), "\\.")[[1]][4], "1x")
  locus_plots_df <- locus_plots_df %>% 
    dplyr::mutate(locus_dir=file.path("www/data",dataset_type,dataset,locus)) %>%
    dplyr::mutate(plot_path=file.path(locus_dir,"plots",paste(locus,"locus_plot",zoom,"png",sep=".")),
                  data_path=file.path(locus_dir,"multi_finemap",paste(locus,"multi_finemap","csv.gz",sep=".")),
                  ld_path=file.path(locus_dir,"LD",gsub("\\.RDS",".csv.gz",basename(path))))
  
  # Arbitrarily use one plot per Locus
  if(!is.null(slice_n)){
    locus_plots_df <- locus_plots_df %>%
      dplyr::group_by(locus, dataset_type, LD_ref, zoom) %>% 
      dplyr::slice(slice_n) 
  }   
  return(data.table::data.table(locus_plots_df))
}


# Only need to run this once beforehand
prepare_files <- function(root="/Volumes/Scizor/Fine_Mapping/Data"){
  ####  Collect plot paths ####
  locus_plots_df <- make_locus_df(root=root) #pattern = "^multiview\\.*"
  
  new_plots <- lapply(1:nrow(locus_plots_df), function(i){
    ROW <- locus_plots_df[i,] 
    dir.create(dirname(ROW$plot_path), showWarnings = F, recursive = T) 
    file.copy(from = ROW$path, 
              to = ROW$plot_path, 
              overwrite = T)
    return(new_path)
  }) %>% unlist() 
  
  
  #### Collect table paths ####
  locus_tables_df <- make_locus_df(root=root,
                                   pattern="*\\.Multi-finemap.tsv.gz")   
  
  locus_tables_df$leadSNP <- parallel::mclapply(1:nrow(locus_tables_df), function(i){
    ROW <- locus_tables_df[i,]
    message(ROW$locus_dir)
    dat <- data.table::fread(ROW$path)
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
    data.table::fwrite(dat, ROW$data_path, row.names = F)  
    leadSNP <- subset(dat, leadSNP)$SNP[1]
    return(leadSNP)
  }, mc.cores = parallel::detectCores()
  ) %>% unlist() 

  
  
  #### Collect LD paths ####
  locus_LD_df <- make_locus_df(root=root,
                               pattern="*UKB_LD.RDS|*1KGphase3_LD.RDS|*1KGphase1_LD.RDS")
   
  ld_paths <- parallel::mclapply(1:nrow(locus_LD_df), function(i){
    ROW <- locus_LD_df[i,]
    locus <- ROW$locus
    message(ROW$locus_dir) 
    # Get lead SNP that's ALSO in LD_matrix 
    if(endsWith(ROW$path, suffix = ".RDS")){
      LD_matrix <- readRDS(ROW$path)
    }
    if(endsWith(ROW$path, suffix = ".RData")){
      load(ROW$path)
    }
    
    # Get the lead SNP for this dataset's locus (ID'ed by locus_dir)
    lead_snp <- subset(locus_tables_df, locus_dir==ROW$locus_dir)$leadSNP
    
    LD_df <-  tryCatch(expr = {  
      LD_df <- data.frame(LD_matrix[, lead_snp])
      LD_df <- cbind(SNP=colnames(LD_matrix), LD_df)
      colnames(LD_df)[2] <- lead_snp
      return(LD_df)
    }, 
    error = function(e){
      LD_df <- data.frame(SNP=colnames(LD_matrix),
                          r=rep(NA,length( colnames(LD_matrix))))
      colnames(LD_df)[2] <-lead_snp
      return(LD_df)
    }, finally = {
      dir.create(dirname(ROW$ld_path), showWarnings = F, recursive = T) 
      data.table::fwrite(LD_df, ROW$ld_path, row.names = F) 
      return(LD_df) 
    }) 
    return(ROW$ld_path)
  }, mc.cores = parallel::detectCores() 
  ) %>% unlist()  
   
 
  
  # Merge paths into one df 
  ## (each table is only totally accurate 
  ##if you searched for files of that particular type)  
  # all_paths <- dplyr::select(locus_plots_df, -c(data_path,ld_path)) %>%
  #   data.table::merge.data.table( dplyr::select(locus_tables_df, c(locus_dir,data_path)), 
  #                                 by="locus_dir",
  #                                 all = T) %>%
  #   data.table::merge.data.table( dplyr::select(locus_LD_df, c(locus_dir,ld_path)), 
  #                                 by="locus_dir",
  #                                 all = T)

  
  all_paths <- data.frame(file_path=list.files(path = "www/data", full.names = T, recursive = T), stringsAsFactors = F) %>%
    dplyr::mutate(study_type=basename(dirname(dirname(dirname(dirname(file_path))))),
                  study=basename(dirname(dirname(dirname(file_path)))),
                  locus_dir=dirname(dirname(file_path)),
                  locus=basename(dirname(dirname(file_path))),
                  file_type=basename(dirname(file_path))) %>%
    dplyr::mutate(zoom=ifelse(file_type=="plots",strsplit(file_path,"\\.")[[1]][3],"1x")) %>%
    data.table::data.table()
  
  saveRDS(all_paths, "www/all_paths.RDS")
}



### Prepare inputs ####
all_paths <- readRDS("www/all_paths.RDS")
all_paths <- subset(all_paths, study!="LRRK2")  # subfolder issue
## Dropdown inputs


                 
                 
## Gather studies
studies <-  unique(all_paths$study)
meta <- label_studies(all_paths = all_paths)
## Gather loci
loci <-  unique(all_paths$locus)
default_study <- if("Nalls23andMe_2019" %in% studies) "Nalls23andMe_2019" else studies[1]
default_locus <- if("BST1" %in% loci) "BST1" else subset(all_paths, study==default_study)$locus[1]
# default_study <- "Ripke_2014"
# default_locus <- "1"
zooms <- unique(all_paths$zoom)[!is.na(unique(all_paths$zoom))]
# input <- list(); input$study <- "Ripke_2014"; input$locus <- "1"; output <- list();

#### UI ####
ui <- fluidPage(
  
  sidebarLayout(position = "left", 
                sidebarPanel(style = "position:fixed; width:25%; height:95vh; z-index:100; overflow-y:auto",
                  div(
                    a(imageOutput(outputId = "echoR_logo", inline = T), href="https://github.com/RajLabMSSM/echolocatoR", target="_blank"),
                    h2("echolocatoR Results Portal"),
                    h4("An interactive database for fine-mapping results generated using",a(em("echolocatoR"), href="https://github.com/RajLabMSSM/echolocatoR", target="_blank"),".") 
                  ), 
                  ## starts interactions
                  hr(),  
                  selectInput(inputId = "study", 
                              label = h3("Study", tipify(a(href="#metdata_header",icon("far fa-question-circle"), style="size: 10px;"), title = "Click for study metadata.", placement = "right", trigger = "hover") ),
                              choices = studies, 
                              selected = default_study,
                  ),  
                  
                  # Create locus options dynamically based on study
                  uiOutput("locus_selection"),  
                  
                  h3("Key"),
                  span(strong("r2", style="display: inline;"),p("Pairwise LD correlation with lead GWAS SNP", style="display: inline;")),br(),
                  span(p("◆",style="color:red; font-size:25px; display: inline"), p("Lead GWAS SNP", style="display: inline;")),br(),
                  span(p("○",style="color:green; font-size:25px; display: inline"), p("Union Credible Set SNP", style="display: inline;")),br(),
                  span(p("○",style="color:goldenrod; font-size:25px; display: inline"), p("Consensus SNP", style="display: inline;")),br(), 
                   
                  h3("Authors"),
                  p(a("Raj Lab", href="www.rajlab.org"),br(),
                    "Dept. of Neuroscience, Dept. of Genetics & Genomic Sciences",br(), 
                    "Icahn School of Medicine at Mount Sinai, New York, NY"),
                  
                  a(imageOutput(outputId = "GH_logo",inline = T), href="https://github.com/RajLabMSSM/Fine_Mapping_Shiny", target="_blank"),
                ),
                
                
                
                #### Plots & Data ####
                mainPanel( 
                  fluidRow(
                    column(width = 12,
                      h3("Interactive plot", tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title = paste("Instructions:","Hover for details.","Click & drag to zoom in. Double click to zoom out."), placement = "right", trigger = "hover") ),    
                      p("GWAS/QTL summary statistics."),
                      # Checkboxes
                      shiny::checkboxInput(inputId = "separate_finemap_methods", 
                                           label = "Separate fine-mapping methods", 
                                           value = FALSE),
                      h4(textOutput("locus_name"), align="center"),
                      h5(textOutput("n_snps"), align="center"),
                      # plotOutput("plot")
                      plotly::plotlyOutput("plotly", 
                                           height = "500px")
                    ) 
                  ), 
                  fluidRow( 
                    column(width = 12, 
                      h3("Static plot"),  
                      # uiOutput("plots", inline = T), 
                      tabsetPanel(id = "tabset",
                        tabPanel(title = "1x",
                          imageOutput("plot_1x",inline = T)
                                  )
                      ),
                      uiOutput("plots",inline = T),
                     
                    )
                  ),
                  br(),
                  fluidRow( 
                    column(width = 12, 
                      h3("Fine-mapping results"), 
                      p("Standardized GWAS/QTL summary statistics and fine-mapping results for the selected locus."),
                      DT::dataTableOutput("results")
                     # tableOutput("results")
                    ) 
                  ),
                  br(),
                  fluidRow( 
                    column(width = 12, 
                           h3("Study metadata", id="metdata_header"), 
                           DT::dataTableOutput("metadata"), 
                    ) 
                  ), 
                  br()
                ),  
  ),
)




#### SERVER ####  
server <- function(input, output, session) {
  #### Logos ####
  output$echoR_logo <- renderImage({ 
    tryCatch(expr = {list(src="./www/icons/echolocatoR_logo-min.png",  width="150px", margin="0px",padding="0px") },
             error=function(e){
               list(src = "./shiny_input/icons/under_construction.png", width="40px")
             })  
  })
  output$GH_logo <- renderImage({ 
    tryCatch(expr = {list(src="./www/icons/github-logo.png",  width="50px", margin="0px",padding="0px") },
             error=function(e){
               list(src = "./shiny_input/icons/under_construction.png", width="40px")
             })  
  })
  
  
  
  ### Dynamically render locus options
  output$locus_selection <- renderUI({ 
    loci <- unique(subset(all_paths, study==input$study)$locus)
    selectInput(inputId = "locus", 
                label = h3("Locus", tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title = "NOTE: Locus names do not necessarily reflect the causal gene(s).", placement = "right", trigger = "hover") ),  
                choices = loci, 
                selected = default_locus,  
    ) 
  })

  #### Gather data + LD ####
  import_data <- function(input){  
    data_path <- subset(all_paths, study==input$study & locus==input$locus & file_type=="multi_finemap")$file_path[1] 
    ld_path <- subset(all_paths, study==input$study & locus==input$locus & file_type=="LD")$file_path[1] 
    finemap_DT <- data.table::fread(data_path) 
    # Import and merged LD info 
    LD_df <- data.table::fread(ld_path) 
    LD_df$r2 <- LD_df[,2]^2 
    finemap_DT <- data.table::merge.data.table(x = finemap_DT,
                                               y = LD_df,
                                               by = "SNP")
    return(finemap_DT)
  }
  
  # Plot   
  
  
  #### Iterate over multiple zoomed static views ####
  ## Reset plots to avoid previous plots appearing
  get_plots_df <- function(input, zoom=NULL){
    plots_sub <- subset(all_paths, study==input$study & locus==input$locus & file_type=="plots")
    if(any(zoom %in% plots_sub$zoom) ){
      plots_sub <- subset(plots_sub, zoom=zoom)
    }
    return(plots_sub)
  }
   
   
  
  
  #### EXAMPLE ITERATIVE PLOTS
  # Source: https://gist.github.com/wch/5436415/
  
  # # Insert the right number of plot output objects into the web page
  # output$plots <- renderUI({
  #   plots_sub <- get_plots_df(input=input)
  #   plot_output_list <- lapply(unique(plots_sub$zoom), function(z) {
  #     plotname <- paste("plot", z, sep="_")
  #     plot_output <- plotOutput(plotname, height = 280, width = 250) 
  #     appendTab(inputId="tabset",
  #               tabPanel(title = verbatimTextOutput(plotname), plot_output)
  #     )
  #   })
  #   # Convert the list to a tagList - this is necessary for the list of items
  #   # to display properly.
  #   do.call(tagList, plot_output_list)
  # })
  
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (z in zooms) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({  
      print(paste("z:",z)) 
      plotname <- paste("plot", z, sep="_")
      print(paste("plotname:",plotname))
      output[[plotname]] <- renderImage({    
        plots_sub <- get_plots_df(input, zoom = z) 
        print(paste("plots_sub:",plots_sub))
        plot_info <- tryCatch(expr = {list(src=plots_sub$file_path, width="100%") },
                              error=function(e){
                                list(src = "./shiny_input/icons/under_construction.png", width="40px")
                              }) 
        return(plot_info) 
      }, deleteFile = F) 
    })
  }
  
  
  
  
  
  
  #### Locus name ####
  output$locus_name <- renderText({paste("Locus:",input$locus)})
  
  #### Interactive plot ####
  output$plotly <- plotly::renderPlotly({
    withProgress(message = 'Loading', value = 0, {
      finemap_DT <- import_data(input)
      finemap_methods <- gsub("\\.PP", "",grep(pattern = "\\.PP", colnames(finemap_DT), value = T ))
      output$n_snps <- renderText({ paste(length(unique(finemap_DT$SNP)),"SNPs") })
      
      n_steps <- 1+1+length(finemap_methods)+1+1
      incProgress(1/n_steps, detail = paste("Importing data...")) 
      
      
      snp_plot <- function(finemap_DT, 
                           y_var="-log10(P)",
                           ylab_prefix=NULL,
                           locus=NULL, 
                           interactive=T,
                           viridis_color=F,
                           ylimits=NULL, 
                           facet_formula=NULL){ 
        if(y_var!="-log10(P)"){
          finemap_DT[[y_var]] <- as.numeric(finemap_DT[[y_var]])
        } 
        snp.labs <- echolocatoR::construct_SNPs_labels(finemap_DT,
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
                     aes_string(x="Mb", y=y_var), 
                     shape=snp.labs$shape, 
                     color=snp.labs$color, 
                     size=snp.labs$size, 
                     stroke=2, 
                     alpha=.7, show.legend = F) + 
          labs(title =locus,# ifelse(is.null(locus), NULL, paste("Locus :", locus)), 
               y=paste(ylab_prefix,y_var)) + 
          scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0,.2))) + 
          theme_bw() +
          theme(axis.title.y = element_text(size=8), 
                plot.title = element_text(hjust = .5)) 
        
        if(!is.null(ylimits)){
          gp <- suppressMessages(gp + ylim(ylimits))
        }
        if(viridis_color){
          gp <- gp + scale_color_viridis_c()
        } else {
          gp <- gp + scale_color_gradient(low = "blue", high = "red", na.value = "gray",
                                          limits=c(0,1), breaks=c(0,.5,1))
        }
        if(interactive){
          return( plotly::ggplotly(gp) )
        } else {return(gp)}
      }
      
      plt_list <- c()
      incProgress(1/n_steps, detail = paste("Creating GWAS plot..."))
      plt_list[["gwas"]] <- snp_plot(finemap_DT, 
                                     y_var="-log10(P)", 
                                     ylab_prefix = "GWAS")
      
      
      #### Separate fine-mapping methods ####
      if(input$separate_finemap_methods){
        for(m in finemap_methods){
          incProgress(1/n_steps, detail = paste("Creating",m,"plot..."))
          print(m)
          y_var <- paste0(m,".PP") 
          plt_list[[m]] <- snp_plot(finemap_DT = finemap_DT[Support>0,], 
                                    y_var=y_var, 
                                    viridis_color =T,
                                    ylimits = c(0,1.1))
        }  
        output$plotly_height <- renderText({"700px"})
      } else {
        output$plotly_height <- renderText({"500px"})
      }  
      pltly <- plotly::subplot(plt_list, 
                               nrows = length(plt_list), 
                               shareY = F, shareX = T,
                               titleY = T)
      incProgress(1/n_steps, detail = paste("Rendering merged plots...")) 
      return(pltly)  
    }) 
  }) # withProgress
  
  
  #### Data #### 
  ## DT ptions
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
  
  # Results DT
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
  
 
  ### Metadata ####
  output$metadata <- DT::renderDT({ 
    DT::datatable(data = meta, 
                  extensions = c("Buttons","Scroller","FixedColumns","FixedHeader","RowGroup"),
                  class='cell-border stripe compact table-hover',
                  options = getOpts(file_name = "study_metadata",
                                    # Can't have rowGroup and FixedColumns at the same time
                                    rowGroup = NULL#list(dataSrc = consensus_index)
                  ),
                  filter='top', 
                  rownames=F, 
                  selection='single')
    })
   
}


shinyApp(ui, server)


