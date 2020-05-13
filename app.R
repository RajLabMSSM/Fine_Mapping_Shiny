
library(shiny)
library(dplyr)
library(DT)
library(data.table)
library(ggplot2)
library(plotly) 


printer <- function(..., v=T){if(v){print(paste(...))}}

construct_SNPs_labels <- function(subset_DT, lead=T, method=T, consensus=T, verbose=F, remove_duplicates=T){ 
  printer("+ PLOT:: Constructing SNP labels...", v=verbose)
  labelSNPs <- data.table::data.table() 
  subset_DT <- data.table::as.data.table(subset_DT)
  
  ## BEFORE fine-mapping  
  if(lead){
    before <- subset( subset_DT %>% arrange(P), leadSNP == T)
    before$type <- "Lead SNP"
    before$color <- "red"
    before$shape <- 18
    before$size <- 3
    labelSNPs <- rbind(labelSNPs, before, fill=T)
  }
  if(method){
    # AFTER fine-mapping
    after = subset(subset_DT, Support>0) 
    if(dim(after)[1]>0){
      after$type <- "Credible Set"  
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
      cons_SNPs$type <- "Consensus SNP"
      cons_SNPs$color <- "darkgoldenrod1"
      cons_SNPs$shape <- 1
      cons_SNPs$size=4
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

make_locus_dict <- function(root="Data/GWAS/Nalls23andMe_2019",
                            pattern="*_ggbio*.png|multi_finemap_plot.png"){ 
  locus_plots_df <- data.frame(path=list.files(path = root, 
                                               pattern = pattern,
                                               full.names = T, recursive = T), stringsAsFactors = F) %>% 
    dplyr::mutate(locus=basename(dirname(dirname(path)))) %>%
    dplyr::arrange(path) %>%
    dplyr::group_by(locus) %>% 
    dplyr::slice(1) 
  locus_plots <- setNames(locus_plots_df$path, locus_plots_df$locus)
  return(locus_plots)
}


# Only need to run this once beforehand
prepare_files <- function(root="~/Desktop/Fine_Mapping/Data/GWAS/Nalls23andMe_2019"){
  # Collect plot paths
  locus_plots <- make_locus_dict(root=root,
                                 pattern="_ggbio.png|multi_finemap_plot.png")  
  dir.create("www/plots", showWarnings = F, recursive = T)
  locus_plots <- lapply(names(locus_plots), function(locus){
    new_path <- file.path("www","plots",paste0(locus,"_plot.png"))
    file.copy(from = locus_plots[[locus]], 
              to = new_path)
    return(new_path)
  }) %>% unlist() %>% `names<-`(names(locus_plots))
  saveRDS(locus_plots, "www/locus_plots.RDS")
  
  # Collect table paths 
  locus_tables <- make_locus_dict(root=root,
                                  pattern="Multi-finemap_results.txt")   
  dir.create("www/data", showWarnings = F, recursive = T)
  locus_tables <- lapply(names(locus_tables), function(locus){
    new_path <- file.path("www","data",paste0(locus,"_data.csv"))
    dat <- data.table::fread(locus_tables[[locus]])
    data.table::fwrite(dat, new_path)
    locus_tables[[locus]] <<- new_path
  }) %>% unlist() %>% `names<-`(names(locus_tables))
  saveRDS(locus_tables, "www/locus_tables.RDS")
}


locus_plots <- readRDS("www/locus_plots.RDS")
# Collect table paths 
locus_tables <- readRDS("www/locus_tables.RDS")

# Dropdown bars
studies <- "Nalls23andMe_2019" #list.dirs("./Data/GWAS/", full.names = F, recursive = F)
loci <-  sort(setNames(names(locus_tables), names(locus_tables) ) )

#### --------- UI ---------- ####
ui <- fluidPage(
  
  sidebarLayout(position = "left", 
                sidebarPanel(style = "position:fixed; width:20%; height:95vh; overflow-y:auto",
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
                  hr(),
                  h4("Interactive plot directions"),
                  shiny::p("Layered results from fine-mapping each locus."), 
                  shiny::p("Hover over each point (SNP) to show more details."),
                  shiny::p("Click and drag to zoom in. Double click to zoom out."),
                  h5("Key:"),
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
                      plotly::plotlyOutput("plotly", height = "600px")
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
      # finemap_DT <- data.table::fread( "./Data/GWAS/Nalls23andMe_2019/WNT3/Multi-finemap/Multi-finemap_results.txt" )
      finemap_DT <- data.table::fread(locus_tables[[input$locus]] )
      finemap_methods <- gsub("\\.PP", "",grep(pattern = "\\.PP", colnames(finemap_DT), value = T ))
      output$n_snps <- renderText({ paste(length(unique(finemap_DT$SNP)),"SNPs") })
      
      n_steps=1+1+length(finemap_methods)+1+1
      incProgress(1/n_steps, detail = paste("Importing data...")) 
      
      finemap_DT$Mb <- finemap_DT$POS/1000000  
      finemap_DT <- find_consensus_SNPs(finemap_DT, verbose = F)
      
      snp_plot <- function(finemap_DT, 
                           y_var="-log10(P)",
                           ylab_prefix=NULL,
                           locus=NULL, 
                           interactive=T,
                           viridis_color=F,
                           ylimits=NULL){
        snp.labs <- construct_SNPs_labels(finemap_DT,
                                          remove_duplicates = F)
        gp <- ggplot(data=finemap_DT, aes(x=Mb, y=eval(parse(text = y_var)), 
                                          color=eval(parse(text = y_var)),
                                          label=CHR,
                                          label1=SNP,
                                          label2=Effect,
                                          label3=P, 
                                          label4=StdErr,
                                          label5=A1, 
                                          label6=A2 ),) + 
          geom_point(alpha=.5, show.legend = F, size=1) +
          geom_point(data = snp.labs, 
                     aes(x=Mb, y=eval(parse(text = y_var))), 
                     shape=snp.labs$shape, 
                     color=snp.labs$color, 
                     size=snp.labs$size - 1, 
                     alpha=.7, show.legend = F) + 
          labs(title =locus,# ifelse(is.null(locus), NULL, paste("Locus :", locus)), 
               y=paste(ylab_prefix,y_var), color=y_var) + 
          scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0,.2)))  +
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
      for(m in finemap_methods){
        incProgress(1/n_steps, detail = paste("Creating",m,"plot..."))
        plt_list[[m]] <- snp_plot(finemap_DT, 
                                  y_var=paste0(m,".PP"), 
                                  viridis_color = F,
                                  ylimits = c(0,1.1))
      }
      
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
  getOpts <- function(file_name="finemapping_results"){
    opts <- list(extensions = c('FixedColumns',"FixedHeader","Buttons"),
                 scrollY = 500, sScrollX="100%", bScrollCollapse=T, pageLength=50,
                 dom = 'frtipB', buttons = list( list(extend = 'csv', filename=file_name),
                                                 list(extend = 'excel', filename=file_name),
                                                 list(extend = 'pdf', filename=file_name),
                                                 list(extend = 'print', filename=file_name),
                                                 'copy'), paging=F,
                 fixedColumns = list(leftColumns = 1))
    return(opts)
  }
  
  # Results table
  output$results <- DT::renderDT({ 
    finemap_DT <- data.table::fread(locus_tables[[input$locus]])
    # createTable(finemap_DT)
    DT::datatable(finemap_DT, 
                  options = getOpts(file_name = paste0(input$locus,"results",sep="_")),
                  filter='top', rownames=F,
                  class='cell-border stripe',
                  selection='single',
                  extensions=c('Buttons','Scroller'))
  })
}


shinyApp(ui, server)


