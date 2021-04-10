
library(shiny)
library(shinyBS)
library(dplyr)
# devtools::install_github('rstudio/DT')
library(DT)
library(data.table)
library(R.utils)
library(ggplot2)
library(plotly) 


# Shiny struggles to find BioC packages...
## https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/3
# library(BiocManager)
# options(repos = BiocManager::repositories())

# To get error logs, run this:
# rsconnect::showLogs()


#### UI ####
ui <- fluidPage(
  
  sidebarLayout(position = "left", 
                #### Header ####
                sidebarPanel(style = "position:fixed; width:25%; height:95vh; z-index:100; overflow-y:auto",
                  div(
                    a(imageOutput(outputId = "echoR_logo", inline = T), href="https://github.com/RajLabMSSM/echolocatoR", target="_blank"),
                    h2("echolocatoR",br(),"Fine-mapping",br(),"Portal"),
                    h4("An interactive database for fine-mapping results generated using",a(em("echolocatoR."), href="https://github.com/RajLabMSSM/echolocatoR", target="_blank")) 
                  ), 
                  #### Study selection ### 
                  uiOutput("study_selection"),   
                  #### Locus selection ####
                  # Create locus options dynamically based on study
                  uiOutput("locus_selection"),  
                  #### LD selection ####
                  # Create LD options dynamically based on study
                  uiOutput("LD_selection"),  
                  
                  
                  h3("SNP Key"),
                  div(style="background-color: rgba(0,0,0,.9); padding: 10px; border-radius: 10px;",
                      span(style="color: white;",strong("r2", style="display: inline;"),p("Pairwise LD correlation with lead GWAS/QTL SNP", style="display: inline; font-size: 12px;")),br(),
                      span(p("◆",style="color:red; font-size:35px; display: inline"), p("Lead GWAS/QTL SNP", style="display: inline; color: white; font-size: 12px;")),br(),
                      span(p("○",style="color:green; font-size:35px; display: inline"), p("Union Credible Set SNP", style="display: inline; color: white; font-size: 12px;")),br(),
                      span(p("○",style="color:goldenrod; font-size:35px; display: inline"), p("Consensus SNP", style="display: inline; color: white; font-size: 12px;")),br(), 
                  ),  
                  br(),
                  span(downloadButton("downloadData_bulk", "Bulk Download"), 
                       tipify(a(icon("far fa-question-circle"), style="size: 10px;"), 
                              title =  "Merge and download the entire database as a single .csv file." )
                  ), 
                   
                  h3("Authors"),
                  p(a("Raj Lab", href="http://www.rajlab.org/", target="_blank"),br(),
                    "Dept. of Neuroscience",br(),
                    "Dept. of Genetics & Genomic Sciences",br(), 
                    "Icahn School of Medicine at Mount Sinai, New York, NY"),
                  
                  h3("Citation"),
                  p(
                    "If you use data/plots from this repository please cite both the orginal GWAS/QTL publication indicated in the",a("metadata",href="#metadata_header"),"and:",br(),
                    em(a("   BM Schilder, J Humphrey, T Raj (2020) echolocatoR: an automated end-to-end statistical and functional genomic fine-mapping pipeline. bioRxiv, 2020.10.22.351221;", href="https://doi.org/10.1101/2020.10.22.351221", target="_blank")) 
                  ),
                  
                  # GH logo and link to repo
                  br(),
                  a(imageOutput(outputId="GH_logo",inline=T), href="https://github.com/RajLabMSSM/Fine_Mapping_Shiny", target="_blank"),
                ), # END LEFT PANEL
                
                
                
                #### Plots & Data ####
                mainPanel( 
                  fluidRow(
                    column(width = 12,
                      h3("Interactive plot", tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title = paste("Instructions:","Hover for details.","Click & drag to zoom in. Double click to zoom out."), placement = "right", trigger = "hover") ),    
                      p("GWAS/QTL summary statistics."),
                      # Options
                      fluidRow(
                        column(width = 4,
                          shiny::checkboxInput(inputId = "separate_finemap_methods", 
                                               label =  h5("Fine-mapped PP", tipify(a(icon("far fa-question-circle"), style="size: 5px;"), title = paste("Show the posterior probabilities (PP) from each fine-mapping tool under the GWAS track."), placement = "right", trigger = "hover") ), 
                                               value = FALSE)
                        ),
                        #### P-value filters ####
                        column(width = 4, 
                               uiOutput("pv_slider", inline = T) 
                               ),
                        column(width = 2,
                               uiOutput("pv_numeric",  inline = T)  
                               )
                      ),
                     
                      # Locus info
                      h4(textOutput("locus_name"), align="center"),
                      h5(textOutput("n_snps"), align="center"), 
                      plotly::plotlyOutput("plotly", 
                                           # height = verbatimTextOutput(textOutput("plotly_height", inline = T))
                                           height = "500px",
                                           )
                    ) 
                  ), 
                  fluidRow( 
                    column(width = 12, 
                      h3("Static plots"),   
                      tabsetPanel(id = "plot_tabset", 
                        tabPanel(title = "1x",
                          imageOutput("plot_1x",inline = T)
                                  ),
                        tabPanel(title = "4x",
                                 imageOutput("plot_4x",inline = T)
                        ),
                        tabPanel(title = "10x",
                                 imageOutput("plot_10x",inline = T)
                        ),
                        uiOutput("plot_tabs"),
                        br(),
                        em("* Note: PLAC-seq interactions for astrocytes are missing from these plots as this data is not available",a("(Nott et. al., 2019).",href="https://science.sciencemag.org/content/366/6469/1134.abstract", target="_blank"))
                      ), 
                    )
                  ),
                  br(),
                  fluidRow( 
                    column(width = 12, 
                      ##### Fine-mapping results #####
                      h3("Fine-mapping results"), 
                      p("Standardized GWAS/QTL summary statistics and fine-mapping results for the selected locus."),
                      h5(a("Column descriptions", href="https://github.com/RajLabMSSM/echolocatoR/tree/dev#multi-finemap-results-files", target="_blank"),
                         tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title = paste("Description of columns in multi-finemap results."), placement = "right", trigger = "hover") ),    
                      span(downloadButton("downloadData", "Download Full Data"), 
                           tipify(a(icon("far fa-question-circle"), style="size: 10px;"), 
                                  title =  "To speed up the app, only Union Credible Set and lead GWAS/QTL SNPs are displayed. Download the full data here." )
                           ), 
                      DT::dataTableOutput("results")
                     # tableOutput("results")
                    ) 
                  ),
                  br(),
                  fluidRow( 
                    column(width = 12, 
                           h3("Study metadata", id="metadata_header"), 
                           DT::dataTableOutput("metadata"), 
                    ) 
                  ), 
                  br()
                ),  
  ),
)





#### SERVER ####  
server <- function(input, output, session) {
  
  #### Support functions ####
  
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
  
  printer <- function(..., v=T){if(v){print(paste(...))}}
  
  
  construct_SNPs_labels <- function(subset_DT,
                                    lead=T,
                                    method=T,
                                    consensus=T,
                                    remove_duplicates=T,
                                    verbose=F){
    printer("+ PLOT:: Constructing SNP labels...", v=verbose)
    labelSNPs <- data.table::data.table()
    subset_DT <- data.table::as.data.table(subset_DT)
    
    ## BEFORE fine-mapping
    if(lead){
      before <- subset( subset_DT %>% dplyr::arrange(P), leadSNP == T)
      before$type <- "Lead"
      before$color <- "red"
      before$shape <- 18
      before$size <- 3
      labelSNPs <- rbind(labelSNPs, before, fill=T)
    }
    if(method){
      # AFTER fine-mapping
      after = subset(subset_DT, Support>0)
      if(dim(after)[1]>0){
        after$type <- "UCS"
        after$color<- "green3"
        after$shape <- 1
        after$size=3
        labelSNPs <- rbind(labelSNPs, after, fill=T)
      }
    }
    if(consensus & "Consensus_SNP" %in% colnames(subset_DT)){
      # Conensus across all fine-mapping tools
      cons_SNPs <- subset(subset_DT, Consensus_SNP==T)
      if(dim(cons_SNPs)[1]>0){
        cons_SNPs$type <- "Consensus"
        cons_SNPs$color <- "darkgoldenrod1"
        cons_SNPs$shape <- 1
        cons_SNPs$size=4
        labelSNPs <- rbind(labelSNPs, cons_SNPs, fill=T)
      }
    }
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
    return(as.data.frame(labelSNPs))
  }
  
  
  
  #### Prepare inputs ####
  all_paths <- readRDS("www/all_paths.RDS")
  # Correct subfolder issue
  all_paths <- subset(all_paths, (study!="LRRK2") &
                        (study_type!="Nalls23andMe_2019"))
 
  # Dropdown inputs  
  ## Gather studies
  studies <-  unique(all_paths$study)
  meta <- label_studies(all_paths = all_paths)
  
  ## Remove ALS until we release the preprint.
  meta <- subset(meta, dataset!="Nicolas_2018_hg38") 
  all_paths <- subset(all_paths, study!="Nicolas_2018_hg38")
  
  ## Gather loci
  loci <-  unique(all_paths$locus)
  ## Gather LD panels
  LD_refs <- unique(all_paths$LD_ref)[!is.na(unique(all_paths$LD_ref))]
  default_study <- if("Nalls23andMe_2019" %in% studies) "Nalls23andMe_2019" else studies[1]
  default_locus <- if("BST1" %in% loci) "BST1" else subset(all_paths, study==default_study)$locus[1]
  default_LD <- if("UKB" %in% LD_refs) "UKB" else subset(all_paths, study==default_study & locus==default_locus)$LD_ref[1]
  # finemap_DT <- readRDS("www/BST1.finemap_DT.RDS")
  
  # default_study <- "Ripke_2014"
  # default_locus <- "1"
  zooms <- unique(all_paths$zoom)[!is.na(unique(all_paths$zoom))]
  # input <- list(); input$study <- "Ripke_2014"; input$locus <- "1"; output <- list();
  
  
  
  
  
  
  #### Logos ####
  output$echoR_logo <- renderImage({ 
    tryCatch(expr = {list(src="./www/icons/echolocatoR_logo-min.png",  width="150px", margin="0px",padding="0px") },
             error=function(e){
               list(src = "./shiny_input/icons/under_construction.png", width="40px")
             })  
  }, deleteFile = F)
  
  output$GH_logo <- renderImage({ 
    tryCatch(expr = {list(src="./www/icons/github-logo.png",  width="50px", margin="0px",padding="0px") },
             error=function(e){
               list(src = "./shiny_input/icons/under_construction.png", width="40px")
             })  
  }, deleteFile = F)
  
 
  output$study_selection <- renderUI({ 
    studies <- unique(all_paths$study)
    selectInput(inputId = "study", 
                label = h3("Study", tipify(a(href="#metadata_header",icon("far fa-question-circle"), style="size: 10px;"), title = "Click for study metadata.", placement = "right", trigger = "hover") ),
                choices = studies, 
                selected = default_study,
    ) 
  })
  
  
  ### Dynamically render locus options
  output$locus_selection <- renderUI({ 
    shiny::validate(
      shiny::need(input$study != "", "Checking study...")
    )
    loci <- unique(subset(all_paths, study==input$study)$locus) 
    selectInput(inputId = "locus", 
                label = h3("Locus", tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title = "NOTE: Locus names do not necessarily reflect the causal gene(s).", placement = "right", trigger = "hover") ),  
                choices = loci, 
                selected = default_locus,  
    ) 
  })
  
  ### Dynamically render locus options
  output$LD_selection <- renderUI({ 
    shiny::validate(
      shiny::need(input$study != "", "Checking study..."),
      shiny::need(input$locus != "", "Checking locus...")
    )
    paths <- subset(all_paths, 
                      study==input$study & 
                      locus==input$locus & 
                      # Only include LD options for which we have data
                      file_type=="multi_finemap") 
    LD_refs <- unique(paths$LD_ref)
    LD_refs <- LD_refs[!is.na(LD_refs)]
    selectInput(inputId = "LD_ref", 
                label = h3("LD panel", tipify(a(icon("far fa-question-circle"), style="size: 10px;"), title =paste( "For some studies, we repeated fine-mapping with a different linkage disequilibrium (LD) panel."), placement = "right", trigger = "hover") ),  
                choices = LD_refs, 
                selected = default_LD,  
    ) 
  })
  
   

  #### Import data + LD ####
  import_data <- function(input){   
    shiny::validate(  
      shiny::need(input$study != "", "Checking study..."),
      shiny::need(input$locus != "", "Checking locus..."),
      shiny::need(input$LD_ref != "", "Checking LD_ref...")
    )
    while(!exists("all_paths")){  
      print(paste("while all_paths",all_paths))
      Sys.sleep(.1)
    }
    data_path <- subset(all_paths, study==input$study & locus==input$locus & LD_ref==input$LD_ref & file_type=="multi_finemap")$file_path[1] 
    ld_path <- subset(all_paths, study==input$study & locus==input$locus & LD_ref==input$LD_ref & file_type=="LD")$file_path[1] 
    print(data_path)
    while(!exists("data_path")){  
      print(paste("while data_path",data_path))
      Sys.sleep(.1)
    }
    
    finemap_DT <- data.table::fread(data_path, nThread = 1)
    # print(nrow(finemap_DT))
    while(!exists("finemap_DT")){
      print(paste("while finemap_DT",nrow(finemap_DT)))
      Sys.sleep(.1)
    }
     
    
    if(!"study" %in% colnames(finemap_DT)) finemap_DT <- cbind(study=input$study, finemap_DT)
    if(!"LD_ref" %in% colnames(finemap_DT)) finemap_DT <- cbind(LD_ref=input$LD_ref, finemap_DT)
    ## Filter p-values 
    shiny::validate(
      shiny::need(input$pval_filter != "", "Checking p-value filter slider..."),
      shiny::need(input$pval_filter_numeric != "", "Checking p-value filter numeric...")
    )
    finemap_DT <- subset(finemap_DT, P<=as.numeric(input$pval_filter))
    
    # Import and merged LD info 
    if(!is.na(ld_path)){
      LD_df <- data.table::fread(ld_path, nThread = 1) 
      LD_df$r2 <- LD_df[,2]^2 
      finemap_DT <- data.table::merge.data.table(x = finemap_DT,
                                                 y = LD_df,
                                                 by = "SNP")  
    }else {
      finemap_DT$r2 <- NA
    }
   
    return(finemap_DT)
  }
  
  # Plot   
  #### Iterate over multiple zoomed static views ####
  ## Reset plots to avoid previous plots appearing
  get_plots_df <- function(input, 
                           zoom=NULL){
    shiny::validate(
      shiny::need(input$study != "", "Checking study..."),
      shiny::need(input$locus != "", "Checking locus..."),
      shiny::need(input$LD_ref != "", "Checking LD_ref...")
    )
    plots_sub <- subset(all_paths, 
                          study==input$study &
                          locus==input$locus & 
                          LD_ref==input$LD_ref &
                          file_type=="plots")
    if(!is.null(zoom)){
      z <- zoom
      plots_sub <- subset(plots_sub, zoom==z)
    }
    return(plots_sub)
  }
   
   
  
  
  #### EXAMPLE ITERATIVE PLOTS
  # Source: https://gist.github.com/wch/5436415/
  
  # # Insert the right number of plot output objects into the web page
  # output$plot_tabs <- renderUI({
  #   plots_sub <- get_plots_df(input=input)
  #   ZOOMS <- unique(plots_sub$zoom)[!is.na( unique(plots_sub$zoom))]
  #   
  #   plot_output_list <- lapply(ZOOMS, function(z) { 
  #     file_path <- get_plots_df(input=input, zoom = z)[1,]$file_path 
  #     # appendTab(inputId="plot_tabset",
  #     #           
  #     #           )
  #     tabPanel(title = z, 
  #              imageOutput(paste("image",z,sep="_"), width="100%"))
  #   
  #   })
  #   # Convert the list to a tagList - this is necessary for the list of items
  #   # to display properly.
  #   do.call(tagList, plot_output_list)
  # })
  
  get_zooms <- function(input){
    shiny::validate(
      shiny::need(input$study != "", "Checking study..."),
      shiny::need(input$locus != "", "Checking locus..."),
      shiny::need(input$LD_ref != "", "Checking LD_ref...")
    )
    zooms <- subset(all_paths,
                    study==input$study &
                      locus==input$locus &
                      LD_ref==input$LD_ref &
                      file_type=="plots")$zoom
    zooms <- zooms[!is.na(zooms)]
    return(zooms)
  }
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
 
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
  for (Z in zooms) { 
    local({
        z <- Z
        output[[paste("plot",z,sep="_")]] <- renderImage({    
          # zooms <- get_zooms(input)
          # plot_name <- paste("plot", z, sep="_")
          # print(paste("plot_name:",plot_name))
          plots_sub <- get_plots_df(input, zoom = z)[1,]
            tryCatch(expr = {list(src=plots_sub$file_path, width="100%") },
                     error=function(e){
                       list(src = "./www/icons/loading.gif", width="40px")
                     })
        }, deleteFile = F); 
    })
  }
    
    
    
   
  
  #### Locus name ####
  output$locus_name <- renderText({
    shiny::validate( 
      shiny::need(input$locus != "", "Checking locus...")
    )
    paste("Locus:",input$locus)
    })
  
  #### Interactive plot ####
  output$plotly <- plotly::renderPlotly({
    tryCatch(expr = {
      
   pltly <- withProgress(message = 'Loading', value = 0, { 
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
        snp.labs <- construct_SNPs_labels(finemap_DT,
                                          remove_duplicates = F)
        color_var <- if(all(is.na(finemap_DT$r2))) NULL else "r2"
        label_A1 <-if("A1" %in% colnames(finemap_DT)) "A1" else NULL
        label_A2 <-if("A2" %in% colnames(finemap_DT)) "A2" else NULL
     
        gp <- ggplot(data=finemap_DT, aes_string(x="Mb", y=y_var, 
                                          color=color_var,
                                          label_CHR="CHR",
                                          label_SNP="SNP",
                                          label_Effect="Effect",
                                          label_P="P", 
                                          label_StdErr="StdErr",
                                          label_A1=label_A1, 
                                          label_A2=label_A2) ) + 
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
      dataset_type <- as.character(subset(meta, dataset==input$study)$dataset_type[1])
      plt_list[[dataset_type]] <- snp_plot(finemap_DT, 
                                     y_var="-log10(P)", 
                                     ylab_prefix = dataset_type)
       
      
      #### Separate fine-mapping methods ####
      shiny::validate(
        shiny::need(is.logical(input$separate_finemap_methods), "Checking separate_finemap_methods...") 
      )
      if(input$separate_finemap_methods){
        for(m in finemap_methods){
          incProgress(1/n_steps, detail = paste("Creating",m,"plot..."))
          print(m)
          y_var <- paste0(m,".PP") 
          plt_list[[m]] <- snp_plot(finemap_DT = finemap_DT, 
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
    }) # end withProgress
   
   
   #### Results DT ####
   # Render within plotly function to avoid importing the data twice
   output$results <- DT::renderDT({  
     # finemap_DT <- import_data(input) 
     shiny::validate( 
       shiny::need(input$locus!="","Checking locus within DT...")
     )
     consensus_index <- grep("Consensus_SNP", colnames(finemap_DT))-1
     UCS_only <- T
     if(UCS_only) finemap_DT <- subset(finemap_DT, Support>0 | leadSNP)
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
   return(pltly)
    },
   error=function(e){  
     # Make a silly little plotly instead of showing an error
     txt <- c("L","O","A","D","I","N","G")
     dat <- data.frame(txt=txt, x=1:length(txt), y=rep(c(1,1),4)[1:length(txt)]) 
     gp_tmp <- ggplot(data = dat, aes(x=x, y=y, label=txt, color=x)) + 
       geom_point(show.legend = F) +
       geom_text(nudge_y = .5) +
       scale_color_gradient(low = "darkslateblue",high = "cyan") +
       theme_void() +
       theme(legend.position='none', 
             axis.line.x = element_blank(),
             axis.line.y = element_blank(), 
             panel.grid.major.y = element_blank()
             ) 
    plotly::ggplotly(gp_tmp, tooltip = NULL, height = 1)
   })## End tryCatch
}) # end  renderPlotly
  
  
   
  
  
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
  
  #### Download full data ####
  # Reactive value for selected dataset ----
  datasetInput <- reactive({
    finemap_DT <- import_data(input)
  })
  
  # datasetInput_bulk <- reactive({
  #   withProgress(message = 'Gathering and merging all results', value = 0, {
  #     dat_paths <- subset(all_paths, file_type=="multi_finemap") 
  #     merged_DT <- lapply(1:nrow(dat_paths), function(i){
  #       ROW <- dat_paths[i,]
  #       dat <- data.table::fread(ROW$file_path, nThread = 1)
  #       dat <- cbind(study=ROW$study, study_type=ROW$study_type, LD_ref=ROW$LD_ref, unique(dat))
  #       incProgress(1/nrow(dat_paths), detail = paste(round(i/nrow(dat_paths),2)*100,"%"))
  #       return(dat)
  #     }) %>% data.table::rbindlist(fill = T)
  #   })
  #   return(merged_DT)
  # })
  datasetInput_bulk <- function(){
    withProgress(message = 'Gathering and merging all results', value = 0, {
      dat_paths <- subset(all_paths, file_type=="multi_finemap") 
      merged_DT <- lapply(1:nrow(dat_paths), function(i){
        ROW <- dat_paths[i,]
        dat <- data.table::fread(ROW$file_path, nThread = 1)
        dat <- cbind(study=ROW$study, study_type=ROW$study_type, LD_ref=ROW$LD_ref, unique(dat))
        incProgress(1/nrow(dat_paths), detail = paste(round(i/nrow(dat_paths),2)*100,"%"))
        return(dat)
      }) %>% data.table::rbindlist(fill = T)
    })
    return(merged_DT)
  }
  
   
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$locus, input$study, input$LD_ref,"multi_finemap","csv",sep = ".")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = F) 
    }
  )
  
  
  output$downloadData_bulk <- downloadHandler(
    filename = function() {
      paste("echolocatoR_Finemapping_Portal","all_data","multi_finemap","csv",sep = ".")
    },
    content = function(file) {
      write.csv(datasetInput_bulk(), file, row.names = F)
    }
  )
  
  #### P-value filters #### 
  output$pv_slider <- renderUI({ 
    sliderInput(inputId = "pval_filter",
                label = h5("p-value filter", tipify(a(icon("far fa-question-circle"), style="size: 5px;"), title = paste("WARNING: Raising the p-value threshold too high can cause the browser to become very slow."), placement = "right", trigger = "hover")),
                min = 0,
                max = 1, 
                value = input$pval_filter_numeric) 
  })
  output$pv_numeric <- renderUI({ 
    numericInput("pval_filter_numeric", 
                 label = NULL, 
                 step = .01,
                 min = 0, 
                 max = 1,
                 value = input$pval_filter)
  })
  
  updateSliderInput(session,"pval_filter", value = .05)
  updateNumericInput(session,"pval_filter_numeric", value = .05 )
  
  
   
}# END SERVER


shinyApp(ui, server)


