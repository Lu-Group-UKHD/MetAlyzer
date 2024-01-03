library(shiny)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)

# setwd('/Users/qianwu/Desktop/RShiny_Biocrates_DataAnalysis')
# metabObj <- MetAlyzer_dataset(file_path = './data/extraction_data.xlsx', silent = T)

ui <- fluidPage(
  # App title
  titlePanel('Biocrates Metabolomics Analysis'),
  
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Preprocessing',
      sidebarLayout(
        sidebarPanel(
          tags$h4('Data uploading', style = 'color:steelblue;font-weight: bold'),
          fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                    placeholder = 'No .xlsx file selected'),
          # Show filtering options only after file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data filtering', style = 'color:steelblue;font-weight: bold'),
                           selectInput('featChoicesFiltering', 'Select feature(s) to remove:',
                                       choices = character(0), multiple = T),
                           selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                       choices = character(0), multiple = T),
                           fluidRow(
                             column(width = 6, actionButton('updateFiltering', 'Filter', width = '100%')),
                             column(width = 6, actionButton('resetFiltering', 'Reset', width = '100%'))
                           ),
                           tags$br(),
                           tags$h4('Log₂(FC) computation', style = 'color:steelblue;font-weight: bold'),
                           fluidRow(
                             column(width = 5, selectInput('smpChoiceGpsLog2FC', 'Compare between:',
                                                           choices = 'Not available', multiple = F)),
                             column(width = 3, selectInput('smpChoicesLog2FC_1', 'Group1',
                                                           choices = character(0), multiple = F),
                                    offset = 1),
                             column(width = 3, selectInput('smpChoicesLog2FC_2', 'Group2',
                                                           choices = character(0), multiple = F))
                           ),
                           checkboxGroupInput('plotLog2FC', 'Visualize:',
                                              choices = c('Volcano plot', 'Scatter plot'))
          )
        ),
        mainPanel(
          fluidRow(
            column(width = 7, verbatimTextOutput('textConcSumm')),
            column(width = 5, verbatimTextOutput('textQuanSumm'))
          ),
          tableOutput('outputId'),
          plotOutput('plotVolcano'),
          plotOutput('plotScatter')
        )
      )
    ),
    tabPanel(
      'Something',
      tags$h4('Could be a section for something else')
    )
  )
)

server <- function(input, output, session) {
  # Create reactive objects for storing up-to-date data
  reactMetabObj <- reactiveValues(metabObj = NULL, ori_metabObj = NULL)
  reactLog2FCTab <- reactiveVal()
  
  # Initialize MetAlyzer SE object
  observeEvent(input$uploadedFile, {
    validUploadedFile <- try(
      metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                    sheet = 1, silent = T),
      silent = T)
    if (!is(validUploadedFile, 'try-error')) {
      reactMetabObj$metabObj <- metabObj
      # Make copy of original data for filtering reset
      reactMetabObj$ori_metabObj <- metabObj
      # Monitor whether uploaded file is valid to show further operations on client side
      output$ifValidUploadedFile <- reactive({
        !is.null(reactMetabObj$metabObj)
      })
      outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
    } else {
      showModal(modalDialog(
        title = 'Uploaded file reading failed...',
        'Is the uploaded .xlsx file generated from MetIDQ™ software?',
        easyClose = T,
        footer = NULL
      ))
      reactMetabObj$metabObj <- NULL
      reactMetabObj$ori_metabObj <- NULL
    }
  })
  
  
  # Prepare choices for feature filtering
  featChoices <- reactive({
    req(reactMetabObj$metabObj)
    list(Class = as.list(unique(rowData(reactMetabObj$metabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactMetabObj$metabObj)))
  })
  # Prepare needed information for sample filtering (i.e., choices, dictionary for
  # choice group searches, and list of choice groups containing duplicated choices)
  # and log2(FC) calculation
  smpChoicePack <- reactive({
    req(reactMetabObj$metabObj)
    smpChoiceList <- as.list(colData(reactMetabObj$metabObj))
    # Uniquify choices in every choice group and convert NA into character so
    # that it can be shown on client side
    for (i in seq_along(smpChoiceList)) {
      choices <- unique(smpChoiceList[[i]])
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        smpChoiceList[[i]] <- c(choices, 'NA')
      } else {
        smpChoiceList[[i]] <- choices
      }
    }
    # Remove choice groups with only one level
    # rmChoiceGps <- which(sapply(smpChoiceList, length) == 1)
    # if (length(rmChoiceGps) != 0) {
    #   smpChoiceList <- smpChoiceList[-rmChoiceGps]
    # }
    # Make identical choices in different choice groups unique by adding suffixes
    allSmpChoices <- unlist(smpChoiceList)
    dupSmpChoices <- allSmpChoices[duplicated(allSmpChoices)]
    # Collect choice groups with duplicated choices for later choice group searches
    # (due to design of MetAlyzer::filterMetaData)
    dupSmpChoiceGps <- c()
    if (length(dupSmpChoices) != 0) {
      for (i in seq_along(smpChoiceList)) {
        choices <- smpChoiceList[[i]]
        if (any(dupSmpChoices %in% choices)) {
          dupIdx <- which(choices %in% dupSmpChoices)
          choices[dupIdx] <- paste0(choices[dupIdx], ' (', names(smpChoiceList)[i], ')')
          smpChoiceList[[i]] <- choices
          dupSmpChoiceGps <- c(dupSmpChoiceGps, names(smpChoiceList)[i])
        }
      }
    }
    
    # Create dictionary for later choice group searches of certain choices (due
    # to design of MetAlyzer::filterMetaData)
    choiceGpSizes <- sapply(smpChoiceList, length)
    smpChoices2Gps <- lapply(seq_along(choiceGpSizes), function(i) {
      rep(names(choiceGpSizes)[i], choiceGpSizes[i])
    })
    smpChoices2Gps <- unlist(smpChoices2Gps)
    names(smpChoices2Gps) <- unlist(smpChoiceList)
    
    return(list(smpChoiceList = smpChoiceList, smpChoices2Gps = smpChoices2Gps,
                dupSmpChoiceGps = dupSmpChoiceGps))
  })
  # Update choices for feature filtering
  observe({
    req(featChoices())
    if ('Metabolism Indicators' %in% unlist(featChoices()) &
        nrow(reactMetabObj$metabObj) == nrow(reactMetabObj$ori_metabObj)) {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
                        selected = 'Metabolism Indicators')
    } else {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
    }
  })
  # Update choices for sample filtering
  observe({
    req(smpChoicePack()$smpChoiceList)
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
  })
  # Do feature filtering
  observeEvent(input$updateFiltering, {
    if (!is.null(input$featChoicesFiltering)) {
      #### Probably add slider for 'min_percent_valid' and options for 'per_group'
      #### to filter unmet features out by threshold of valid values
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = input$featChoicesFiltering)
    } else {
      #### Some features will still be filtered out
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = NULL)
    }
  })
  # Do sample filtering
  observeEvent(input$updateFiltering, {
    req(input$smpChoicesFiltering)
    for (selectedChoice in input$smpChoicesFiltering) {
      selectedChoiceGp <- smpChoicePack()$smpChoices2Gps[selectedChoice]
      # Prepare up-to-date choices of selected choice group to avoid error that
      # metadata of same sample is selected multiple times at once
      choices <- unique(colData(reactMetabObj$metabObj)[[selectedChoiceGp]])
      # Convert NA into character so that it can be searched
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        choices <- c(choices, 'NA')
      }
      if (!selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
        if (selectedChoice %in% choices) {
          if (selectedChoice == 'NA') {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !.data[[selectedChoiceGp]] %in% selectedChoice)
          }
        }
      } else {
        # Prepare manually added suffix if selected choice group is in collected
        # choice groups with duplicates (see smpChoicePack() section)
        regexSuffix <- paste0(' \\(', selectedChoiceGp, '\\)')
        modSelectedChoice <- stringr::str_remove(selectedChoice, regexSuffix)
        if (modSelectedChoice %in% choices) {
          if (modSelectedChoice == 'NA') {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                                                !.data[[selectedChoiceGp]] %in% modSelectedChoice)
          }
        }
      }
    }
  })
  # Reset feature and sample filtering
  observeEvent(input$resetFiltering, {
    reactMetabObj$metabObj <- reactMetabObj$ori_metabObj
  })
  
  
  # Update sample choice groups for log2(FC) calculation
  observe({
    req(smpChoicePack()$smpChoiceList)
    smpChoiceList <- smpChoicePack()$smpChoiceList
    # Remove choice groups whose level sizes are one or equal to sample size
    choiceGpSizes <- sapply(smpChoiceList, length)
    rmChoiceGps <- which(choiceGpSizes == 1 | choiceGpSizes == ncol(reactMetabObj$metabObj))
    if (length(rmChoiceGps) != 0) {
      smpChoiceList <- smpChoiceList[-rmChoiceGps]
    }
    smpChoiceGps <- names(smpChoiceList)
    if (length(smpChoiceGps) != 0) {
      updateSelectInput(session, 'smpChoiceGpsLog2FC', choices = smpChoiceGps)
    } else {
      updateSelectInput(session, 'smpChoiceGpsLog2FC', choices = 'Not available')
    }
  })
  # Update sample choices according to selected choice group for log2(FC) calculation
  observe({
    req(smpChoicePack()$smpChoiceList)
    selectedChoiceGp <- input$smpChoiceGpsLog2FC
    smpChoices <- smpChoicePack()$smpChoiceList[[selectedChoiceGp]]
    if (selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
      regexSuffix <- paste0(' \\(', selectedChoiceGp, '\\)')
      smpChoices <- stringr::str_remove(smpChoices, regexSuffix)
    }
    # Turn NULL into empty vector, or choices will not be updated, i.e., previous
    # choices are shown
    if (length(smpChoices) == 0) {
      smpChoices <- character(0)
    }
    updateSelectInput(session, 'smpChoicesLog2FC_1', choices = smpChoices)
    updateSelectInput(session, 'smpChoicesLog2FC_2', choices = smpChoices, selected = smpChoices[2])
  })
  # Compute log2(FC) of features
  observe({
    req(input$plotLog2FC)
    selectedChoiceGp <- input$smpChoiceGpsLog2FC
    selectedChoices <- c(input$smpChoicesLog2FC_1, input$smpChoicesLog2FC_2)
    if (!'NA' %in% selectedChoices) {
      metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                            .data[[selectedChoiceGp]] %in% selectedChoices)
    } else {
      metabObj <- MetAlyzer::filterMetaData(reactMetabObj$metabObj,
                                            is.na(.data[[selectedChoiceGp]]) |
                                              .data[[selectedChoiceGp]] %in% selectedChoices)
    }
    log2FCTab <- MetAlyzer::calculate_log2FC(metabObj, .data[[selectedChoiceGp]],
                                             perc_of_min = 0.2, impute_NA = T)
    reactLog2FCTab(log2FCTab)
  })
  
  
  # Show statistics of data concentration values and quantification statuses 
  output$textConcSumm <- renderPrint({
    req(reactMetabObj$metabObj)
    MetAlyzer::summarizeConcValues(reactMetabObj$metabObj)
  })
  output$textQuanSumm <- renderPrint({
    req(reactMetabObj$metabObj)
    MetAlyzer::summarizeQuantData(reactMetabObj$metabObj)
  })
  
  output$outputId <- renderTable({
    dim(reactMetabObj$metabObj)
  })
}

shinyApp(ui = ui, server = server)
