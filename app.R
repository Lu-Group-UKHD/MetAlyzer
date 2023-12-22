library(shiny)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)

# setwd('/Users/qianwu/Desktop/RShiny_Biocrates_DataAnalysis')
# metabObj <- MetAlyzer_dataset(file_path = './data/extraction_data.xlsx', silent = T)

ui <- fluidPage(
  # App title
  titlePanel('Biocrates Data Analysis'),
  
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Preprocessing',
      sidebarLayout(
        sidebarPanel(
          fileInput('uploadedFile', 'Upload xlsx file:', multiple = F, accept = '.xlsx'),
          # Show filtering options only after a file is uploaded
          conditionalPanel(condition = "output.if_uploadedFile",
                           selectInput('featFiltering', 'Select feature(s) to remove:',
                                       choices = character(0), multiple = T),
                           selectInput('smpFiltering', 'Select sample(s) to remove:',
                                       choices = character(0), multiple = T),
                           actionButton('updateFiltering', 'Filter', width = '49%'),
                           actionButton('resetFiltering', 'Reset', width = '49%'))
        ),
        mainPanel(
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
  
  # Monitor whether file is uploaded to show further operations on client side
  output$if_uploadedFile <- reactive({
    !is.null(input$uploadedFile)
  })
  outputOptions(output, 'if_uploadedFile', suspendWhenHidden = F)
  
  # Initialize MetAlyzer SE object
  observeEvent(input$uploadedFile, {
    reactObjs$metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                            sheet = 1, silent = T)
    # Make copy of original data for filtering reset
    reactObjs$ori_metabObj <- reactObjs$metabObj
  })
  
  
  # Prepare choices for feature filtering
  featChoices <- reactive({
    list(Class = as.list(unique(rowData(reactObjs$metabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactObjs$metabObj)))
    
  })
  # Prepare choices for sample filtering
  smpChoices <- reactive({
    tmp_smpChoices <- as.list(colData(reactObjs$metabObj))
    # Uniquify choices in every choice group and convert NA into character so
    # that it can be shown on client side
    for (i in seq_len(length(tmp_smpChoices))) {
      choices <- unique(tmp_smpChoices[[i]])
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        tmp_smpChoices[[i]] <- c(choices, 'NA')
      } else {
        tmp_smpChoices[[i]] <- choices
      }
    }
    # Remove choice groups with only one level
    # rm_choiceGps <- which(sapply(tmp_smpChoices, length) == 1)
    # if (length(rm_choiceGps) != 0) {
    #   tmp_smpChoices <- tmp_smpChoices[-rm_choiceGps]
    # }
    # Make identical choices in different choice groups unique by adding suffixes
    all_smpChoices <- unlist(tmp_smpChoices)
    dup_smpChoices <- all_smpChoices[duplicated(all_smpChoices)]
    # Collect choice groups with duplicated choices for later choice group searches
    # (due to design of MetAlyzer::filterMetaData)
    dup_smpChoiceGps <- c()
    if (length(dup_smpChoices) != 0) {
      for (i in seq_len(length(tmp_smpChoices))) {
        choices <- tmp_smpChoices[[i]]
        if (any(dup_smpChoices %in% choices)) {
          dupIdx <- which(choices %in% dup_smpChoices)
          choices[dupIdx] <- paste0(choices[dupIdx], ' (', names(tmp_smpChoices)[i], ')')
          tmp_smpChoices[[i]] <- choices
          dup_smpChoiceGps <- c(dup_smpChoiceGps, names(tmp_smpChoices)[i])
        }
      }
    }
    # Make choices lists so that sole choice in certain choice group can be shown
    tmp_smpChoices <- lapply(tmp_smpChoices, function(choices) {
      as.list(choices)
    })
    
    # Create dictionary for later choice group searches of certain choices (due
    # to design of MetAlyzer::filterMetaData)
    choiceGpSizes <- sapply(tmp_smpChoices, length)
    tmp_smpChoices2Gps <- lapply(seq_len(length(choiceGpSizes)), function(i) {
      rep(names(choiceGpSizes)[i], choiceGpSizes[i])
    })
    tmp_smpChoices2Gps <- unlist(tmp_smpChoices2Gps)
    names(tmp_smpChoices2Gps) <- unlist(tmp_smpChoices)
    
    return(list(smpChoices = tmp_smpChoices, smpChoices2Gps = tmp_smpChoices2Gps,
                dup_smpChoiceGps = dup_smpChoiceGps))
  })
  # Update choices for feature and sample filtering
  observe({
    req(input$uploadedFile)
    # For feature filtering
    if ('Metabolism Indicators' %in% unlist(featChoices()) &
        nrow(reactObjs$metabObj) == nrow(reactObjs$ori_metabObj)) {
      updateSelectInput(session, 'featFiltering', choices = featChoices(),
                        selected = 'Metabolism Indicators')
    } else {
      updateSelectInput(session, 'featFiltering', choices = featChoices())
    }
    # For sample filtering
    updateSelectInput(session, 'smpFiltering', choices = smpChoices()$smpChoices)
  })
  # Do feature and sample filtering
  observeEvent(input$updateFiltering, {
    # For feature filtering
    #### Probably add slider for 'min_percent_valid' and options for 'per_group'
    #### to filter unmet features out by threshold of valid values
    reactObjs$metabObj <- MetAlyzer::filterMetabolites(reactObjs$metabObj,
                                                       drop_metabolites = input$featFiltering)
    # For sample filtering
    for (selectChoice in input$smpFiltering) {
      selectChoiceGp <- smpChoices()$smpChoices2Gps[selectChoice]
      if (!selectChoiceGp %in% smpChoices()$dup_smpChoiceGps) {
        if (selectChoice == 'NA') {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          !is.na(.data[[selectChoiceGp]]))
        } else {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          .data[[selectChoiceGp]] != selectChoice)
        }
      } else {
        # Prepare manually added suffix if selected choice group is in collected
        # choice groups with duplicates (see smpChoices() section)
        regexSuffix <- paste0(' \\(', selectChoiceGp, '\\)')
        tmp_selectChoice <- stringr::str_remove(selectChoice, regexSuffix)
        if (tmp_selectChoice == 'NA') {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          !is.na(.data[[selectChoiceGp]]))
        } else {
          reactObjs$metabObj <- MetAlyzer::filterMetaData(reactObjs$metabObj,
                                                          .data[[selectChoiceGp]] != tmp_selectChoice)
        }
      }
    }
  })
  # Reset feature and sample filtering
  observeEvent(input$resetFiltering, {
    reactObjs$metabObj <- reactObjs$ori_metabObj
  })
  
  
  output$outputId <- renderTable({
    dim(reactObjs$metabObj)
  })
}

shinyApp(ui = ui, server = server)
