library(shiny)
library(shinycssloaders)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)
source("utils.R")

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
          tags$h4('Data uploading', style = 'color:steelblue;font-weight:bold'),
          fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                    placeholder = 'No .xlsx file selected'),
          # Show filtering options only after file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data filtering', style = 'color:steelblue;font-weight:bold'),
                           selectInput('featChoicesFiltering', 'Select metabolite(s) to remove:',
                                       choices = character(0), multiple = T),
                           selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                       choices = character(0), multiple = T),
                           fluidRow(
                             column(width = 6, actionButton('updateFiltering', 'Filter', width = '100%')),
                             column(width = 6, actionButton('resetFiltering', 'Reset', width = '100%'))
                           ),
          )
        ),
        mainPanel(
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data distribution', style = 'color:Black;font-weight:bold'),
                           plotlyOutput('plotDatDist')  %>%
                             shinycssloaders::withSpinner(color="#56070C"),
                           tags$h4('Sample metadata', style = 'color:Black;font-weight:bold'),
                           DT::dataTableOutput('tblSmpMetadat')  %>%
                             shinycssloaders::withSpinner(color="#56070C"),
                           tags$h4('Quantification status', style = 'color:Black;font-weight:bold'),
                           plotlyOutput('plotQuanStatus')  %>%
                             shinycssloaders::withSpinner(color="#56070C"),
          )
        )
      )
    ),
    tabPanel(
      'Visualisation',
      sidebarLayout(
        sidebarPanel(
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Log₂(FC) Calculation', style = 'color:steelblue;font-weight:bold'),
                           fluidRow(
                             column(width = 5, selectInput('smpChoiceGpsLog2FC', 'Compare between:',
                                                           choices = 'Not available', multiple = F)),
                             column(width = 3, selectInput('smpChoicesLog2FC_1', 'Group1',
                                                           choices = character(0), multiple = F),
                                    offset = 1),
                             column(width = 3, selectInput('smpChoicesLog2FC_2', 'Group2',
                                                           choices = character(0), multiple = F))
                           ),
                           fluidRow(
                             column(width = 5, actionButton('computeLog2FC', 'Compute', width = '100%'))
                           ),
                           tags$br(),
                           tags$h4('Log₂(FC) Visualization', style = 'color:steelblue;font-weight:bold'),
                           fluidRow(
                             column(width = 4, checkboxInput('plotVulcanoLog2FC', 'Vulcano Plot',
                                                             value = FALSE, width = '100%')),
                             column(width = 4, checkboxInput('plotScatterLog2FC', 'Scatter Plot',
                                                             value = FALSE, width = '100%')),
                             column(width = 4, checkboxInput('plotNetworkLog2FC', 'Network Plot',
                                                             value = FALSE, width = '100%'))
                           ),
                           fluidRow(
                             style = "display: flex; align-items: center;",
                             column(width = 7, selectInput('metabChoicesVulcano',
                                                           'Select metabolite(s) to highlight:',
                                                           choices = character(0), multiple = T)),
                             column(width = 5, checkboxInput('highlightVulcano', 'Highlight',
                                                             value = FALSE, width = '100%'))
                           ),
          )
        ),
        mainPanel(
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           # Use conditionalPanel to make selected plot always shown at top
                           conditionalPanel(condition = "input.plotVulcanoLog2FC",
                                            plotlyOutput('plotVolcano') %>%
                                              shinycssloaders::withSpinner(color="#56070C"),
                           ),
                           #### Can spinner be moved?
                           conditionalPanel(condition = "input.plotScatterLog2FC",
                                            fluidRow(
                                              column(width = 9, plotlyOutput('plotScatter') %>%
                                                       shinycssloaders::withSpinner(color="#56070C")),
                                              column(width = 3, plotOutput('plotScatterLegend'))
                                            ),
                           ),
                           conditionalPanel(condition = "input.plotNetworkLog2FC",
                                            plotOutput('plotNetwork') %>%
                                              shinycssloaders::withSpinner(color="#56070C")
                           ),
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Create reactive objects for storing up-to-date data
  reactMetabObj <- reactiveValues(metabObj = NULL, ori_metabObj = NULL)
  reactLog2FCTbl <- reactiveVal()
  reactVulcanoHighlight <- reactiveVal()
  
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
  
  # Retrieve abundance data and sample metadata for showing data overviews
  datOverviewTbls <- reactive({
    req(reactMetabObj$metabObj)
    metabAggreDat <- metadata(reactMetabObj$metabObj)$aggregated_data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    metabSmpMetadat <- colData(reactMetabObj$metabObj) %>%
      tibble::as_tibble(rownames = 'ID') %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    # Use original column names whose spaces are not replaced with '.'
    colnames(metabSmpMetadat) <- c('ID', colnames(colData(reactMetabObj$metabObj)))
    # Prepare ID levels for displaying samples in order
    idLevels <- rownames(colData(reactMetabObj$metabObj))
    metabAggreDat <- dplyr::left_join(metabAggreDat, metabSmpMetadat, by = 'ID') %>%
      dplyr::mutate(ID = factor(ID, levels = paste0('Smp', idLevels)))
    
    return(list(metabAggreDat = metabAggreDat, metabSmpMetadat = metabSmpMetadat))
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
                                                             drop_metabolites = input$featChoicesFiltering,
                                                             drop_NA_concentration = NULL)
    } else {
      #### Some features will still be filtered out
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = NULL,
                                                             drop_NA_concentration = NULL)
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
    # Turn NULL (when selectedChoiceGp is 'Not available') into empty vector, or
    # choices will not be updated, i.e., previous choices are shown
    if (length(smpChoices) == 0) {
      smpChoices <- character(0)
    }
    updateSelectInput(session, 'smpChoicesLog2FC_1', choices = smpChoices)
    updateSelectInput(session, 'smpChoicesLog2FC_2', choices = smpChoices, selected = smpChoices[2])
  })
  # Compute log2(FC)
  observeEvent(input$computeLog2FC, {
    if (input$smpChoicesLog2FC_1 != input$smpChoicesLog2FC_2) {
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
      # Use do.call to prepare arguments to be passed due to design of function:
      # deparse(substitute(categorical))
      log2FCTbl <- do.call(MetAlyzer::calculate_log2FC, list(metalyzer_se = metabObj,
                                                             categorical = as.symbol(selectedChoiceGp),
                                                             impute_perc_of_min = 0.2, impute_NA = T))
      reactLog2FCTbl(log2FCTbl)
    } else {
      showModal(modalDialog(
        title = 'Log₂(FC) computation failed...',
        'Please select two different sample groups.',
        easyClose = T,
        footer = NULL
      ))
    }
  })

  
  # Update choices for metabolite highlighting
  observe({
    req(featChoices())
    updateSelectInput(session, 'metabChoicesVulcano', choices = featChoices())
  })
  # Create new column for highlighting metabolites
  #### Probably better change 'highlight_metabolites' to 'highlight'
  observeEvent(input$highlightVulcano, {
    req(reactLog2FCTbl())
    if (!is.null(input$metabChoicesVulcano)) {
      selectedChoices <- input$metabChoicesVulcano
      highlightTbl <- reactLog2FCTbl()
      
      highlightTbl$highlight_metabolites <- "Other Metabolites"
      highlightTbl$highlight_metabolites[highlightTbl$Metabolite %in% selectedChoices] <- "Highlighted Metabolite(s)"
      highlightTbl$highlight_metabolites <- as.factor(highlightTbl$highlight_metabolites)
      
      reactVulcanoHighlight(highlightTbl)
    }
  })
  observeEvent(input$metabChoicesVulcano, {
    req(reactLog2FCTbl())
    if (input$highlightVulcano) {
      selectedChoices <- input$metabChoicesVulcano
      highlightTbl <- reactLog2FCTbl()
      
      highlightTbl$highlight_metabolites <- "Other Metabolites"
      highlightTbl$highlight_metabolites[highlightTbl$Metabolite %in% selectedChoices] <- "Highlighted Metabolite(s)"
      highlightTbl$highlight_metabolites <- as.factor(highlightTbl$highlight_metabolites)
      
      reactVulcanoHighlight(highlightTbl)
    }
  })
  
  
  # Set plot theme
  # th <- theme_bw(base_size = 15) +
  #   theme(axis.title = element_text(face = 'bold'),
  #         axis.text = element_text(face = 'bold'),
  #         axis.ticks = element_line(linewidth = 0.8),
  #         legend.text = element_text(size = 15))
  
  # Show data overviews
  # Data distribution
  output$plotDatDist <- renderPlotly({
    req(datOverviewTbls()$metabAggreDat)
    metabAggreDat <- datOverviewTbls()$metabAggreDat
    ggplotly(
      ggplot(metabAggreDat, aes(x=ID, y=Concentration)) +
        geom_boxplot() +
        scale_y_log10() +
        labs(x = 'Sample', y = 'Metabolite abundance') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  # Sample metadata
  output$tblSmpMetadat <- DT::renderDataTable({
    req(datOverviewTbls()$metabSmpMetadat)
    metabSmpMetadat <- datOverviewTbls()$metabSmpMetadat
    DT::datatable(metabSmpMetadat, rownames = F, filter = list(position = 'top', clear = T, plain = F),
                  selection = list(mode = 'single', target = 'row'), style = 'bootstrap')
  })
  # Quantification status
  output$plotQuanStatus <- renderPlotly({
    req(datOverviewTbls()$metabAggreDat)
    metabAggreDat <- datOverviewTbls()$metabAggreDat %>%
      dplyr::group_by(ID, Status) %>%
      dplyr::summarise(Count = dplyr::n())
    # Prepare colors for quantification statuses
    #### Choose better colors and add color for NA quantification status
    quanStatus2Color <- c(Valid = '#1f78b4', LOQ = '#e31a1c', LOD = '#ff7f00', Invalid = '#33a02c')
    quanStatusCols <- quanStatus2Color[names(quanStatus2Color) %in% unique(metabAggreDat$Status)]
    ggplotly(
      ggplot(metabAggreDat, aes(x = ID, y = Count, fill = Status)) +
        geom_col(position = "stack") +
        scale_fill_manual('Status', values = quanStatusCols) +
        labs(x = "Sample", y = "Count") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  
  # Visualize log2(FC)
  # Vulcano plot
  #### Set cutoffs for log2(FC) and p-value as parameters
  #### Highlight also metabolic classes
  #### Use white background by theme_bw()
  output$plotVolcano <- renderPlotly({
    req(reactLog2FCTbl())
    if (input$plotVulcanoLog2FC) {
      if (input$highlightVulcano) {
        req(reactVulcanoHighlight())
        plots <- plotly_log2FC(reactVulcanoHighlight(), vulcano = T, scatter = F)
        plots$HighlightedVulcanoPlot
      } else {
        plots <- plotly_log2FC(reactLog2FCTbl(), vulcano = T, scatter = F)
        plots$VulcanoPlot
      }
    }
  })
  # Scatter plot
  output$plotScatter <- renderPlotly({
    req(reactLog2FCTbl())
    if (input$plotScatterLog2FC) {
      plot <- plotly_log2FC(reactLog2FCTbl(), vulcano = F, scatter = T)
      plot$Scatterplot$Plot
    }
  })
  output$plotScatterLegend <- renderPlot({
    req(reactLog2FCTbl())
    if (input$plotScatterLog2FC) {
      plot <- plotly_log2FC(reactLog2FCTbl(), vulcano = F, scatter=T)
      plot$Scatterplot$Legend
    }
  })
  # Network plot
  output$plotNetwork <- renderPlot({
    req(reactLog2FCTbl())
    if (input$plotNetworkLog2FC) {
      MetAlyzer::plot_network(reactLog2FCTbl())
    }
  })
}

shinyApp(ui = ui, server = server)
