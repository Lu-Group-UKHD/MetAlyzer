library(shiny)
library(shinyBS)
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
  sidebarLayout(
    sidebarPanel(
      tags$h4('Data uploading', style = 'color:steelblue;font-weight:bold'),
      fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                placeholder = 'No .xlsx file selected'),
      # Show filtering options only after file is uploaded
      conditionalPanel(condition = "output.ifValidUploadedFile",
                       tags$h4('Data filtering', style = 'color:steelblue;font-weight:bold'),
                       bsCollapse(
                         open = 'Sample filtering', multiple = T,
                         bsCollapsePanel('Sample filtering', style = 'info',
                                         selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                                     choices = character(0), multiple = T)),
                         bsCollapsePanel(
                           'Metabolite filtering', style = 'info',
                           fluidRow(
                             column(width = 8, sliderInput('featCompleteCutoffFiltering',
                                                           'Select % of observed values each feature should have:',
                                                           min = 0, max = 100, value = 80))
                           ),
                           fluidRow(
                             column(width = 8, sliderInput('featValidCutoffFiltering',
                                                           'Select % of valid values each feature should have:',
                                                           min = 0, max = 100, value = 67)),
                             column(width = 3, offset = 1,
                                    checkboxGroupInput('featValidStatusFiltering', 'Validity',
                                                       choices = c('Valid', 'LOQ', 'LOD', 'Invalid'),
                                                       selected = c('Valid', 'LOQ')))
                           )
                           # selectInput('featChoicesFiltering', 'Select metabolite(s) to remove:',
                           #             choices = character(0), multiple = T)
                         )
                       ),
                       fluidRow(
                         column(width = 6, actionButton('updateFiltering', 'Filter', width = '100%')),
                         column(width = 6, actionButton('resetFiltering', 'Reset', width = '100%'))
                       ),
                       tags$br(),
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
      tabsetPanel(
        type = 'tabs',
        tabPanel(
          'Preprocessing',
        
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           bsCollapse(
                             open = 'Data distribution', multiple = T,
                             bsCollapsePanel('Sample metadata', style = 'primary',
                                             DT::dataTableOutput('tblSmpMetadat') %>%
                                               withSpinner(color="#56070C")),
                             bsCollapsePanel('Data distribution', style = 'primary',
                                             plotlyOutput('plotDatDist') %>%
                                               withSpinner(color="#56070C")),
                             bsCollapsePanel('Data completeness', style = 'primary',
                                             bsCollapsePanel('Hint', style = 'success',
                                                             textOutput('summDatComplete', container = strong)),
                                             plotlyOutput('plotDatComplete') %>%
                                               withSpinner(color="#56070C")),
                             bsCollapsePanel('Quantification status', style = 'primary',
                                             #### Adjust tab size
                                             bsCollapsePanel('Hint', style = 'success',
                                                             textOutput('summQuanStatus', container = strong)),
                                             plotlyOutput('plotQuanStatus') %>%
                                               withSpinner(color="#56070C"))
                           ),
          )
        ),
        tabPanel(
          'Visualisations',
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           # Use conditionalPanel to make selected plot always shown at top
                           conditionalPanel(condition = "input.plotVulcanoLog2FC",
                                            plotlyOutput('plotVolcano') %>%
                                              withSpinner(color="#56070C"),
                           ),
                           #### Can spinner be moved?
                           conditionalPanel(condition = "input.plotScatterLog2FC",
                                            fluidRow(
                                              column(width = 9, style = "z-index:2;", plotlyOutput('plotScatter') %>%
                                                       shinycssloaders::withSpinner(color="#56070C")),
                                              column(width = 3, style = "margin-left: -175px; z-index:1;",
                                                     imageOutput('plotScatterLegend'))
                                            ),
                           ),
                           conditionalPanel(condition = "input.plotNetworkLog2FC",
                                            plotlyOutput('plotNetwork') %>%
                                              shinycssloaders::withSpinner(color="#56070C")
                            )
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
  
  # Retrieve abundance data and sample metadata and compute feature completeness
  # level and quantification status validity level for showing data overviews
  datOverviewPack <- reactive({
    req(reactMetabObj$metabObj)
    metabAggreTbl <- metadata(reactMetabObj$metabObj)$aggregated_data %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    smpMetadatTbl <- colData(reactMetabObj$metabObj) %>%
      tibble::as_tibble(rownames = 'ID') %>%
      dplyr::mutate(ID = paste0('Smp', ID))
    # Use original column names whose spaces are not replaced with '.'
    colnames(smpMetadatTbl) <- c('ID', colnames(colData(reactMetabObj$metabObj)))
    # Prepare ID levels for displaying samples in order
    idLevels <- rownames(colData(reactMetabObj$metabObj))
    metabAggreTbl <- dplyr::left_join(metabAggreTbl, smpMetadatTbl, by = 'ID') %>%
      dplyr::mutate(ID = factor(ID, levels = paste0('Smp', idLevels)))
    
    # Compute completeness level of each feature for doing filtering
    featCompleteLvTbl <- dplyr::select(metabAggreTbl, Metabolite, Concentration) %>%
      dplyr::mutate(Obs = dplyr::case_when(!Concentration %in% c(0, NA) ~ 1),
                    Miss = dplyr::case_when(Concentration %in% c(0, NA) ~ 1)) %>%
      dplyr::group_by(Metabolite) %>%
      dplyr::summarise(ObsCount = sum(Obs, na.rm = T), MissCount = sum(Miss, na.rm = T)) %>%
      dplyr::mutate(TotalCount = ncol(reactMetabObj$metabObj),
                    CompleteRatio = ObsCount / TotalCount)
    # Compute status counts of each feature for doing filtering
    featStatusCountTbl <- dplyr::select(metabAggreTbl, Metabolite, Status) %>%
      dplyr::mutate(Valid = dplyr::case_when(Status %in% 'Valid' ~ 1),
                    LOQ = dplyr::case_when(Status %in% 'LOQ' ~ 1),
                    LOD = dplyr::case_when(Status %in% 'LOD' ~ 1),
                    Invalid = dplyr::case_when(Status %in% 'Invalid' ~ 1)) %>%
      dplyr::group_by(Metabolite) %>%
      dplyr::summarise(ValidCount = sum(Valid, na.rm = T), LOQCount = sum(LOQ, na.rm = T),
                       LODCount = sum(LOD, na.rm = T), InvalidCount = sum(Invalid, na.rm = T)) %>%
      dplyr::mutate(TotalCount = ncol(reactMetabObj$metabObj))
    
    return(list(metabAggreTbl = metabAggreTbl, smpMetadatTbl = smpMetadatTbl,
                featCompleteLvTbl = featCompleteLvTbl, featStatusCountTbl = featStatusCountTbl))
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
  # observe({
  #   req(featChoices())
  #   if ('Metabolism Indicators' %in% unlist(featChoices()) &
  #       nrow(reactMetabObj$metabObj) == nrow(reactMetabObj$ori_metabObj)) {
  #     updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
  #                       selected = 'Metabolism Indicators')
  #   } else {
  #     updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
  #   }
  # })
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
  #### Conduct sample filtering prior to feature filtering, results in more reliable filtering
  observeEvent(input$updateFiltering, {
    req(datOverviewPack())
    # Collect features to remove based on missingness
    featCompleteLvTbl <- datOverviewPack()$featCompleteLvTbl
    featCompleteLevels <- featCompleteLvTbl$CompleteRatio
    featCompleteCutoff <- input$featCompleteCutoffFiltering / 100
    rmMissFeats <- featCompleteLvTbl$Metabolite[featCompleteLevels < featCompleteCutoff]
    # Collect features to remove based on validity
    featStatusCountTbl <- datOverviewPack()$featStatusCountTbl
    featValidCutoff <- input$featValidCutoffFiltering / 100
    featValidStatus <- input$featValidStatusFiltering
    featValidCounts <- rep(0, nrow(featStatusCountTbl))
    for (status in featValidStatus) {
      featValidCounts <- featValidCounts + featStatusCountTbl[[paste0(status, 'Count')]]
    }
    featValidLevels <- featValidCounts / featStatusCountTbl[['TotalCount']]
    rmInvalidFeats <- featStatusCountTbl$Metabolite[featValidLevels < featValidCutoff]
    # Summarize features to remove
    rmFeats <- unique(c(rmMissFeats, rmInvalidFeats))
    if (length(rmFeats) != 0) { #!is.null(input$featChoicesFiltering)
      #### Check if result is same as using 'min_percent_valid' and 'valid_status'
      reactMetabObj$metabObj <- MetAlyzer::filterMetabolites(reactMetabObj$metabObj,
                                                             drop_metabolites = rmFeats,
                                                             drop_NA_concentration = NULL,
                                                             min_percent_valid = NULL,
                                                             valid_status = c('Valid', 'LOQ'),
                                                             per_group = NULL)
    } else {
      #### Some features will still be filtered out due to 'drop_NA_concentration'.
      #### Set it to NULL, instead of FALSE
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
      metabObj <- MetAlyzer::calculate_log2FC(metalyzer_se = metabObj,
                                                             categorical = selectedChoiceGp,
                                                             impute_perc_of_min = 0.2, impute_NA = T)
      reactLog2FCTbl(log2FC(metabObj))
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
      highlightTbl$highlight_metabolites[highlightTbl$Class %in% selectedChoices] <- "Highlighted Metabolite(s)"
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
      highlightTbl$highlight_metabolites[highlightTbl$Class %in% selectedChoices] <- "Highlighted Metabolite(s)"
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
    req(datOverviewPack()$metabAggreTbl)
    metabAggreTbl <- datOverviewPack()$metabAggreTbl
    ggplotly(
      ggplot(metabAggreTbl, aes(x=ID, y=Concentration)) +
        geom_boxplot() +
        scale_y_log10() +
        labs(x = 'Sample', y = 'Metabolite abundance') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  
  # Sample metadata
  output$tblSmpMetadat <- DT::renderDataTable({
    req(datOverviewPack()$smpMetadatTbl)
    smpMetadatTbl <- datOverviewPack()$smpMetadatTbl
    DT::datatable(smpMetadatTbl, rownames = F, filter = list(position = 'top', clear = T, plain = F),
                  selection = list(mode = 'single', target = 'row'), style = 'bootstrap')
  })
  
  # Data completeness
  output$summDatComplete <- renderText({
    req(datOverviewPack()$featCompleteLvTbl)
    featCompleteLevels <- datOverviewPack()$featCompleteLvTbl$CompleteRatio
    paste(sum(featCompleteLevels < 0.8), 'out of', nrow(reactMetabObj$metabObj),
          'metabolites fail to fulfil 80% rule, which is recommended filtering out.',
          'Besides, is there any sample with high level of missingness?')
  })
  output$plotDatComplete <- renderPlotly({
    req(datOverviewPack()$metabAggreTbl)
    smpCompleteCountTbl <- datOverviewPack()$metabAggreTbl %>%
      dplyr::mutate(Completeness = dplyr::case_when(!Concentration %in% c(0, NA) ~ 'Observed',
                                                    Concentration %in% c(0, NA) ~ 'Missing')) %>%
      dplyr::group_by(ID, Completeness) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    ggplotly(
      ggplot(smpCompleteCountTbl, aes(x=ID, y=Count, fill=Completeness)) +
        geom_col(position = 'stack') +
        scale_fill_manual(values = c('grey', 'black')) +
        labs(x = 'Sample') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  
  # Quantification status
  output$summQuanStatus <- renderText({
    req(datOverviewPack()$featStatusCountTbl)
    featStatusValid <- datOverviewPack()$featStatusCountTbl %>%
      dplyr::mutate(TotalValidCount = ValidCount + LOQCount,
                    ValidRatio = TotalValidCount / TotalCount)
    paste(sum(featStatusValid$ValidRatio < 0.67), 'out of', nrow(reactMetabObj$metabObj),
          'metabolites have less than 67% measurements with valid status (Valid, LOQ),',
          'which is recommended filtering out. Besides, is there any sample containing few valid values?')
  })
  output$plotQuanStatus <- renderPlotly({
    req(datOverviewPack()$metabAggreTbl)
    smpStatusCountTbl <- datOverviewPack()$metabAggreTbl %>%
      dplyr::group_by(ID, Status) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    # Prepare colors for quantification statuses
    #### Add color for NA quantification status
    status2Color <- c(Valid = '#33a02c', LOQ = '#1f78b4', LOD = '#ff7f00', Invalid = '#e31a1c')
    statusCols <- status2Color[names(status2Color) %in% unique(smpStatusCountTbl$Status)]
    ggplotly(
      ggplot(smpStatusCountTbl, aes(x = ID, y = Count, fill = Status)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = statusCols) +
        labs(x = "Sample") +
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
        plotly_vulcano(reactVulcanoHighlight())
      } else {
        plotly_vulcano(reactLog2FCTbl())
      }
    }
  })
  
  # Scatter plot
  output$plotScatter <- renderPlotly({
    req(reactLog2FCTbl())
    if (input$plotScatterLog2FC) {
      plot <- plotly_scatter(reactLog2FCTbl())
      plot$Plot
    }
  })
  output$plotScatterLegend <- renderImage({
    req(reactLog2FCTbl()) 
    if (input$plotScatterLog2FC) {
      plot <- plotly_scatter(reactLog2FCTbl())
      legend <- plot$Legend
      # A temp file to save the output.
      # This file will be removed later by renderImage

      outfile <- tempfile(fileext='.svg')
      ggsave(file=outfile, plot=legend)

      # Return a list containing the filename
      list(src = normalizePath(outfile),
          contentType = 'image/svg+xml',
          width = 571.675,
          height = 400,
          alt = "My svg Histogram")
    }
  })
  
  # Network plot
  output$plotNetwork <- renderPlotly({
    req(reactLog2FCTbl())
    if (input$plotNetworkLog2FC) {
      plotly_network(reactLog2FCTbl())
    }
  })
}

shinyApp(ui = ui, server = server)
