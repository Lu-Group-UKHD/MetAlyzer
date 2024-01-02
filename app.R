library(shiny)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)

# setwd('/Users/qianwu/Desktop/RShiny_Biocrates_DataAnalysis')
# metabObj <- MetAlyzer_dataset(file_path = './data/extraction_data.xlsx', silent = T)

ui <- fluidPage(
  # App title
  titlePanel('Biocrates Metabolomics Data Analysis'),
  
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Preprocessing',
      sidebarLayout(
        sidebarPanel(
          fileInput('uploadedFile', 'Upload xlsx file:', multiple = F, accept = '.xlsx'),
          # Show filtering options only after a file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           selectInput('featFiltering', 'Select feature(s) to remove:',
                                       choices = character(0), multiple = T),
                           selectInput('smpFiltering', 'Select sample(s) to remove:',
                                       choices = character(0), multiple = T),
                           actionButton('updateFiltering', 'Filter', width = '49%'),
                           actionButton('resetFiltering', 'Reset', width = '49%'))
        ),
        mainPanel(
          fluidRow(
            column(width = 7, verbatimTextOutput('textConcSumm')),
            column(width = 5, verbatimTextOutput('textQuanSumm'))
          ),
          tableOutput('outputId')
        )
      )
    ),
    tabPanel(
      'Visualization',
      tags$h4('Could be a section for data visualization')
    )
  )
)

server <- function(input, output, session) {
  # Create reactive object for storing uploaded file
  reactObjs <- reactiveValues(metabObj = NULL, ori_metabObj = NULL)
  
  # Initialize MetAlyzer SE object
  observeEvent(input$uploadedFile, {
    validUploadedFile <- try(
      metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                    sheet = 1, silent = T),
      silent = T)
    if (!is(validUploadedFile, 'try-error')) {
      reactObjs$metabObj <- metabObj
      # Make copy of original data for filtering reset
      reactObjs$ori_metabObj <- metabObj
      # Monitor whether uploaded file is valid to show further operations on client side
      output$ifValidUploadedFile <- reactive({
        !is.null(reactObjs$metabObj)
      })
      outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
    } else {
      showModal(modalDialog(
        title = 'Uploaded file reading failed...',
        'Is the uploaded xlsx file generated from the MetIDQ™ software?',
        easyClose = T,
        footer = NULL
      ))
      reactObjs$metabObj <- NULL
      reactObjs$ori_metabObj <- NULL
    }
  })
  
  
  # Prepare choices for feature filtering
  featChoices <- reactive({
    req(reactObjs$metabObj)
    list(Class = as.list(unique(rowData(reactObjs$metabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactObjs$metabObj)))
  })
  # Prepare needed information for sample filtering (i.e., choices, dictionary for
  # choice group searches, and list of choice groups containing duplicated choices)
  smpFilterList <- reactive({
    req(reactObjs$metabObj)
    smpChoices <- as.list(colData(reactObjs$metabObj))
    # Uniquify choices in every choice group and convert NA into character so
    # that it can be shown on client side
    for (i in seq_len(length(smpChoices))) {
      choices <- unique(smpChoices[[i]])
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        smpChoices[[i]] <- c(choices, 'NA')
      } else {
        smpChoices[[i]] <- choices
      }
    }
    # Remove choice groups with only one level
    # rmChoiceGps <- which(sapply(smpChoices, length) == 1)
    # if (length(rmChoiceGps) != 0) {
    #   smpChoices <- smpChoices[-rmChoiceGps]
    # }
    # Make identical choices in different choice groups unique by adding suffixes
    allSmpChoices <- unlist(smpChoices)
    dupSmpChoices <- allSmpChoices[duplicated(allSmpChoices)]
    # Collect choice groups with duplicated choices for later choice group searches
    # (due to design of MetAlyzer::filterMetaData)
    dupSmpChoiceGps <- c()
    if (length(dupSmpChoices) != 0) {
      for (i in seq_len(length(smpChoices))) {
        choices <- smpChoices[[i]]
        if (any(dupSmpChoices %in% choices)) {
          dupIdx <- which(choices %in% dupSmpChoices)
          choices[dupIdx] <- paste0(choices[dupIdx], ' (', names(smpChoices)[i], ')')
          smpChoices[[i]] <- choices
          dupSmpChoiceGps <- c(dupSmpChoiceGps, names(smpChoices)[i])
        }
      }
    }
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoices <- lapply(smpChoices, function(choices) {
      as.list(choices)
    })
    
    # Create dictionary for later choice group searches of certain choices (due
    # to design of MetAlyzer::filterMetaData)
    choiceGpSizes <- sapply(smpChoices, length)
    smpChoices2Gps <- lapply(seq_len(length(choiceGpSizes)), function(i) {
      rep(names(choiceGpSizes)[i], choiceGpSizes[i])
    })
    smpChoices2Gps <- unlist(smpChoices2Gps)
    names(smpChoices2Gps) <- unlist(smpChoices)
    
    return(list(smpChoices = smpChoices, smpChoices2Gps = smpChoices2Gps,
                dupSmpChoiceGps = dupSmpChoiceGps))
  })
  # Update choices for feature filtering
  observe({
    req(featChoices())
    if ('Metabolism Indicators' %in% unlist(featChoices()) &
        nrow(reactObjs$metabObj) == nrow(reactObjs$ori_metabObj)) {
      updateSelectInput(session, 'featFiltering', choices = featChoices(),
                        selected = 'Metabolism Indicators')
    } else {
      updateSelectInput(session, 'featFiltering', choices = featChoices())
    }
  })
  # Update choices for sample filtering
  observe({
    req(smpFilterList()$smpChoices)
    updateSelectInput(session, 'smpFiltering', choices = smpFilterList()$smpChoices)
  })
  # Do feature filtering
  observeEvent(input$updateFiltering, {
    if (!is.null(input$featFiltering)) {
      #### Probably add slider for 'min_percent_valid' and options for 'per_group'
      #### to filter unmet features out by threshold of valid values
      reactObjs$metabObj <- MetAlyzer::filterMetabolites(reactObjs$metabObj,
                                                         drop_metabolites = input$featFiltering)
    } else {
      reactObjs$metabObj <- MetAlyzer::filterMetabolites(reactObjs$metabObj,
                                                         drop_metabolites = NULL)
    }
  })
  # Do sample filtering
  observeEvent(input$updateFiltering, {
    req(input$smpFiltering)
    for (selectChoice in input$smpFiltering) {
      selectChoiceGp <- smpFilterList()$smpChoices2Gps[selectChoice]
      if (!selectChoiceGp %in% smpFilterList()$dupSmpChoiceGps) {
        if (selectChoice == 'NA') {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          !is.na(.data[[selectChoiceGp]]))
        } else {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          .data[[selectChoiceGp]] != selectChoice)
        }
      } else {
        # Prepare manually added suffix if selected choice group is in collected
        # choice groups with duplicates (see smpFilterList() section)
        regexSuffix <- paste0(' \\(', selectChoiceGp, '\\)')
        modSelectChoice <- stringr::str_remove(selectChoice, regexSuffix)
        if (modSelectChoice == 'NA') {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          !is.na(.data[[selectChoiceGp]]))
        } else {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          .data[[selectChoiceGp]] != modSelectChoice)
        }
      }
    }
  })
  # Reset feature and sample filtering
  observeEvent(input$resetFiltering, {
    reactObjs$metabObj <- reactObjs$ori_metabObj
  })
  
  
  # Show statistics of data concentration values and quantification statuses 
  output$textConcSumm <- renderPrint({
    req(reactObjs$metabObj)
    MetAlyzer::summarizeConcValues(reactObjs$metabObj)
  })
  output$textQuanSumm <- renderPrint({
    req(reactObjs$metabObj)
    MetAlyzer::summarizeQuantData(reactObjs$metabObj)
  })
  
  output$outputId <- renderTable({
    dim(reactObjs$metabObj)
  })
}

shinyApp(ui = ui, server = server)
