library(shiny)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)

# setwd('/Users/qianwu/Desktop/RShiny_Biocrates_DataAnalysis')
# metabObj <- MetAlyzer_dataset(file_path = './data/extraction_data.xlsx', silent = T)

ui <- shiny::fluidPage(
  # App title
  shiny::titlePanel('Biocrates Data Analysis'),
  
  shiny::tabsetPanel(
    type = 'tabs',
    shiny::tabPanel(
      'Preprocessing',
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::fileInput('uploadedFile', 'Upload xlsx file:',
                           multiple = F, accept = '.xlsx'),
          shiny::selectInput('featFiltering', 'Select feature(s) to remove:',
                             choices = character(0), multiple = T),
          shiny::selectInput('smpFiltering', 'Select sample(s) to remove:',
                             choices = character(0), multiple = T),
          shiny::actionButton('updateFiltering', 'Filter', width = '49%'),
          shiny::actionButton('resetFiltering', 'Reset', width = '49%')
        ),
        shiny::mainPanel(
          shiny::tableOutput('outputId'),
          shiny::textOutput('testing')
        )
      )
    ),
    shiny::tabPanel(
      'Visualization',
      tags$h4('Could be a section for data visualization')
    )
  )
)

server <- function(input, output, session) {
  # Create reactive variables
  reactObjs <- shiny::reactiveValues(metabObj = NULL, ori_metabObj = NULL)
  
  # Initialize MetAlyzer SE object
  shiny::observeEvent(input$uploadedFile, {
    reactObjs$metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                            sheet = 1, silent = T)
    # Make a copy of original data for filtering reset
    reactObjs$ori_metabObj <- reactObjs$metabObj
  })
  
  # Update choices for feature and sample filtering
  observe({
    # Stop proceeding if file has not been uploaded. Current calculation (reactive
    # expression) and any callers on it will be stopped silently (avoid error when
    # app is just initiated)
    shiny::req(input$uploadedFile)
    # Feature filtering choices
    featChoices <- list(Class = unique(rowData(reactObjs$metabObj)$metabolic_classes),
                        Metabolite = rownames(reactObjs$metabObj))
    shiny::updateSelectInput(session, 'featFiltering', choices = featChoices)
    # Sample filtering choices
    smpChoices <- lapply(as.list(colData(reactObjs$metabObj)), function(metaVar) {
      metaVar <- unique(metaVar)
      if (all(!is.na(metaVar))) {
        metaVar
      } else {
        # Convert NA into a character so that it can be shown on client side
        #### Not only are levels of metadata variables needed but also metadata
        #### variables. Same levels of different variables (e.g., 'NA') have to
        #### be shown. Try having two select boxes using 'conditionalPanel'?
        metaVar <- metaVar[-which(is.na(metaVar))]
        c(metaVar, 'NA')
      }
    })
    # Remove options with only one level
    smpChoices <- smpChoices[-which(sapply(smpChoices, length) == 1)]
    shiny::updateSelectInput(session, 'smpFiltering', choices = smpChoices)
  })
  
  # Do feature and sample filtering
  shiny::observeEvent(input$updateFiltering, {
    #### There are many more filtering criteria in function to explore
    if (!is.null(reactObjs$metabObj)) {
      # Feature filtering
      reactObjs$metabObj <- MetAlyzer::filterMetabolites(reactObjs$metabObj,
                                                         drop_metabolites = input$featFiltering)
      # Sample filtering
      # if (input$smpFiltering == 'NA') {
      #   reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj, is.na(input$ABC))
      # }
    }
  })
  # Reset feature and sample filtering
  shiny::observeEvent(input$resetFiltering, {
    reactObjs$metabObj <- reactObjs$ori_metabObj
  })
  
  output$outputId <- shiny::renderTable({
    dim(reactObjs$metabObj)
  })
  
  output$testing <- shiny::renderText({
    input$smpFiltering
  })
}

shinyApp(ui = ui, server = server)
