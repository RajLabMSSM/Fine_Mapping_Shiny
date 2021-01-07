library(shiny)

ui <- shinyUI(fluidPage(
  
  # Important! : JavaScript functionality to add the Tabs
  tags$head(tags$script(HTML("
    /* In coherence with the original Shiny way, tab names are created with random numbers. 
       To avoid duplicate IDs, we collect all generated IDs.  */
    var hrefCollection = [];

    Shiny.addCustomMessageHandler('addTabToTabset', function(message){
      var hrefCodes = [];
      /* Getting the right tabsetPanel */
      var tabsetTarget = document.getElementById(message.tabsetName);

      /* Iterating through all Panel elements */
      for(var i = 0; i < message.titles.length; i++){
        /* Creating 6-digit tab ID and check, whether it was already assigned. */
        do {
          hrefCodes[i] = Math.floor(Math.random()*100000);
        } 
        while(hrefCollection.indexOf(hrefCodes[i]) != -1);
        hrefCollection = hrefCollection.concat(hrefCodes[i]);

        /* Creating node in the navigation bar */
        var navNode = document.createElement('li');
        var linkNode = document.createElement('a');

        linkNode.appendChild(document.createTextNode(message.titles[i]));
        linkNode.setAttribute('data-toggle', 'tab');
        linkNode.setAttribute('data-value', message.titles[i]);
        linkNode.setAttribute('href', '#tab-' + hrefCodes[i]);

        navNode.appendChild(linkNode);
        tabsetTarget.appendChild(navNode);
      };

      /* Move the tabs content to where they are normally stored. Using timeout, because
         it can take some 20-50 millis until the elements are created. */ 
      setTimeout(function(){
        var creationPool = document.getElementById('creationPool').childNodes;
        var tabContainerTarget = document.getElementsByClassName('tab-content')[0];

        /* Again iterate through all Panels. */
        for(var i = 0; i < creationPool.length; i++){
          var tabContent = creationPool[i];
          tabContent.setAttribute('id', 'tab-' + hrefCodes[i]);

          tabContainerTarget.appendChild(tabContent);
        };
      }, 100);
    });
    "))),
  # End Important
  
  tabsetPanel(id = "mainTabset", 
              tabPanel("InitialPanel1", "Some Text here to show this is InitialPanel1", 
                       actionButton("goCreate", "Go create a new Tab!"),
                       textOutput("creationInfo")
              ),
              tabPanel("InitialPanel2", "Some Text here to show this is InitialPanel2 and not some other Panel")
  ),
  
  # Important! : 'Freshly baked' tabs first enter here.
  uiOutput("creationPool", style = "display: none;")
  # End Important
  
))

server <- function(input, output, session){
  
  # Important! : creationPool should be hidden to avoid elements flashing before they are moved.
  #              But hidden elements are ignored by shiny, unless this option below is set.
  output$creationPool <- renderUI({})
  outputOptions(output, "creationPool", suspendWhenHidden = FALSE)
  # End Important
  
  # Important! : This is the make-easy wrapper for adding new tabPanels.
  addTabToTabset <- function(Panels, tabsetName){
    titles <- lapply(Panels, function(Panel){return(Panel$attribs$title)})
    Panels <- lapply(Panels, function(Panel){Panel$attribs$title <- NULL; return(Panel)})
    
    output$creationPool <- renderUI({Panels})
    session$sendCustomMessage(type = "addTabToTabset", message = list(titles = titles, tabsetName = tabsetName))
  }
  # End Important 
  
  # From here: Just for demonstration
  output$creationInfo <- renderText({
    paste0("The next tab will be named NewTab", input$goCreate + 1)
  })
  
  observeEvent(input$goCreate, {
    nr <- input$goCreate
    
    newTabPanels <- list(
      tabPanel(paste0("NewTab", nr), 
               actionButton(paste0("Button", nr), "Some new button!"), 
               textOutput(paste0("Text", nr))
      ), 
      tabPanel(paste0("AlsoNewTab", nr), sliderInput(paste0("Slider", nr), label = NULL, min = 0, max = 1, value = 1))
    )
    
    output[[paste0("Text", nr)]] <- renderText({
      if(input[[paste0("Button", nr)]] == 0){
        "Try pushing this button!"
      } else {
        paste("Button number", nr , "works!")
      }
    })
    
    addTabToTabset(newTabPanels, "mainTabset")
  })
}

shinyApp(ui, server)