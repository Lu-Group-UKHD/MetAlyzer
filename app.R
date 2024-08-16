library(shiny)
library(shinyBS)
library(shinyWidgets)
library(shinycssloaders)
library(MetAlyzer)
library(SummarizedExperiment)
library(tidyverse)
library(limma)
source("utils.R")
library(bslib)
library(htmlwidgets)
library(svglite)
library(writexl)
library(vsn)

ui <- fluidPage(
  titlePanel('Biocrates Metabolomics Analysis'),
  tabsetPanel(
    type = 'tabs',
    tabPanel(
      'Overview',
      sidebarLayout(
        sidebarPanel(
          tags$h4('Data Uploading', style = 'color:steelblue;font-weight:bold'),
          uiOutput('updateFileInput'),
          div(textOutput('textFileInput'), style = 'color:IndianRed;font-weight:bold;font-size:110%'),
          checkboxInput('exampleFile',
                        HTML('Explore app with <b>example dataset</b>: <a href = "https://doi.org/10.3389/fmolb.2022.961448">[Gegner et al. 2022]</a>'), 
                        value = FALSE),
          bsTooltip('exampleFile', 'Discover the full range of functions with this diverse tissue dataset.'),
          # Show data processing options only after file is uploaded
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           tags$h4('Data Processing', style = 'color:steelblue;font-weight:bold'),
                           bsCollapse(
                             open = c('Sample filtering', 'Filtering Log'), multiple = T,
                             bsCollapsePanel('Sample filtering', style = 'info',
                                             selectInput('smpChoicesFiltering', 'Select sample(s) to remove:',
                                                         choices = character(0), multiple = T)),
                             bsCollapsePanel(
                               'Metabolite filtering', style = 'info',
                               selectInput('featChoicesFiltering', 'Select metabolite(s) to remove:',
                                           choices = character(0), multiple = T),
                               fluidRow(
                                 column(width = 8, sliderInput('featCompleteCutoffFiltering',
                                                               'Select % of observed values each metabolite should have:',
                                                               min = 0, max = 100, value = 80))
                               ),
                               fluidRow(
                                 column(width = 8, sliderInput('featValidCutoffFiltering',
                                                               'Select % of valid values each metabolite should have:',
                                                               min = 0, max = 100, value = 50)),
                                 column(width = 3, offset = 1,
                                        checkboxGroupInput('featValidStatusFiltering', 'Validity',
                                                           choices = c('Valid', 'LOQ', 'LOD', 'Invalid'),
                                                           selected = c('Valid', 'LOQ'))),
                                 bsTooltip('featCompleteCutoffFiltering',
                                           'Metabolites with observed values below this cutoff are removed.'),
                                 bsTooltip('featValidCutoffFiltering',
                                           'Metabolites with valid values below this cutoff are removed.'),
                                 bsTooltip('featValidStatusFiltering',
                                           'The selected is considered valid for filtering.')
                               )
                             ),
                             bsCollapsePanel(
                               'Imputation and Normalization', style = 'info',
                               materialSwitch('imputation', 'Half-minimum (HM) imputation',
                                              value = T, status = 'primary', right = T),
                               #### Median normalization has to be fixed!
                               selectInput('normalization', 'Select normalization method to use:',
                                           choices = c('None', 'Median log₂ normalization',
                                                       'Total ion count (TIC) log₂ normalization'),
                                           selected = 'Median log₂ normalization', multiple = F)
                               # bsTooltip('imputation', paste('Missing values are replaced with half of the minimum of,
                               #                               observed values in each metabolite.')),
                               # bsTooltip('normalization', 'Median scaling is conducted followed by log2 transformation.')
                             )
                           ),
                           fluidRow(
                             column(width = 6, actionButton('updateProcessing', 'Process', width = '100%')),
                             column(width = 6, actionButton('revertProcessing', 'Revert/Default', width = '100%')),
                             bsTooltip('revertProcessing', 'The processed data and specified parameters revert to the origins.')
                           )
                           #### Hide following functionalities until work is done
                           # tags$br(),
                           # tags$h4('Abundance Matrix Download', style = 'color:steelblue;font-weight:bold'),
                           # fluidRow(
                           #   style = "display: flex; align-items: flex-end;",
                           #   column(width = 7, selectInput('metaChoicesRawAbunExport', 'Select sample metadata to include:',
                           #                                 choices = character(0), multiple = T)),
                           #   column(width = 5, downloadButton('downloadRawAbunExport', 'Download',
                           #                                    style = 'width:100%; margin-bottom: 15px')),
                           #   bsTooltip('downloadRawAbunExport', 'This output can directly be used for MetaboAnalyst.')
                           # ),
                           # tags$br(),
                           # tags$h4('Metabolite Identifier Download', style = 'color:steelblue;font-weight:bold'),
                           # fluidRow(
                           #   style = "display: flex; align-items: flex-end;",
                           #   column(width = 7, uiOutput('loadFeatIdChoicesExport')),
                           #   column(width = 5, downloadButton('downloadFeatIdsExport', 'Download',
                           #                                    style = 'width:100%; margin-bottom: 15px')),
                           #   #### Biocrates names in mapping table and example dataset do not greatly overlap (164 out of 630)?
                           #   #### Create mapping table using materials Gernot provided?
                           #   bsTooltip('featIdChoicesExport', 'The identifiers were generated by MetaboAnalyst using HMDB IDs.',
                           #             placement = 'top')
                           # )
          )
        ),
        mainPanel(
          # JavaScript: Double equals (==) is loose equality with type coercion and
          # triple equals (===) is strict equality without type coercion
          conditionalPanel(condition = "output.ifValidUploadedFile === undefined",
                           HTML('<br>
                                <h3>This tool is here to help you analyze and explore your Biocrates data.</h3>
                                <br>
                                <p>Please ensure your Excel file contains the cell <b>"Class"</b> and the column 
                                <b>"Sample Type"</b> with correct indentation, as these are used as anchor cells!</p>
                                <br>
                                <p>Additional information is available by hovering over all action elements.</p>')
          ),
          conditionalPanel(condition = "output.ifValidUploadedFile",
                           bsCollapse(
                             open = c('Data distribution', 'Data completeness'), multiple = T,
                             bsCollapsePanel('Sample metadata', style = 'primary',
                                             DT::dataTableOutput('tblSmpMetadat') %>%
                                               withSpinner(color="#56070C")),
                             bsCollapsePanel('Data distribution', style = 'primary',
                                             uiOutput('updateGpColsDatDist'),
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
                           )
          )
        )
      )
    ), # TabPanel 1 End
    tabPanel(
      'Log₂(FC)',
      sidebarLayout(
        conditionalPanel(condition = "output.ifValidUploadedFile",
                         sidebarPanel(
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
                           tags$h4('Vulcano Plot', style = 'color:steelblue;font-weight:bold'), #Log₂(FC) Visualization
                           #### Highlighting in scatter plot is to be fixed 
                           fluidRow(
                             style = "display: flex; align-items: flex-end;",
                             column(width = 6, selectInput('metabChoicesVulcano',
                                                           'Select metabolite(s) to highlight:',
                                                           choices = character(0), multiple = T)),
                             column(width = 6, materialSwitch('highlightVulcano', 'Highlight',
                                                              value = F, status = 'primary'))
                           ),
                           tags$h4('Select cutoffs:', style = 'font-weight:bold;font-size:14px'),
                           fluidRow(
                             column(width = 6, sliderInput('plotVolcanoLog2FCCutoff', 'Log₂(FC)',
                                                           min = 0, max = 10, value = 1, step = 0.1, ticks = F)),
                             column(width = 6, selectInput('plotVolcanoPValCutoff', 'P-value',
                                                           choices = c('0.0001', '0.001', '0.01', '0.05', '0.1'), #to avoid scientific notation
                                                           multiple = F, selected = 0.05))
                           )
                         )
        ),
        mainPanel(
          conditionalPanel(condition = "output.ifValidUploadedFile & input.computeLog2FC == 0",
                           div(textOutput('textLog2FC'), style = 'color:IndianRed;font-weight:bold;font-size:110%')),
          conditionalPanel(condition = "input.computeLog2FC",
                           tags$h4(strong('Vulcano plot'), style = "margin-top:1rem;"),
                           plotlyOutput('plotVolcano') %>%
                             withSpinner(color="#56070C"),
                           #### Provide option of downloading static plot
                           fluidRow(style="display:flex; justify-content:right; margin-top:1rem;",
                                    column(width = 2, downloadButton("downloadVulcanoPlot", "Download vulcano plot"))),
                           tags$br(),
                           tags$h4(strong('Scatter plot'), style = "margin-top:1rem;"),
                           fluidRow(
                             column(width = 9, style = "z-index:2;", plotlyOutput('plotScatter') %>%
                                      shinycssloaders::withSpinner(color="#56070C")),
                             column(width = 3, style = "margin-left: -175px; z-index:1;",
                                    imageOutput('plotScatterLegend'))
                           ),
                           fluidRow(style="display:flex; justify-content:right; margin-top:1rem; margin-bottom:1rem;",
                                    column(width = 2, downloadButton("downloadScatterPlot", "Download scatter plot")))
          )
        )
      )
    ), # TabPanel 2 End
    tabPanel(
      #### Add title to scale, i.e., log2(FC)?
      'Network',
      # Lower network plot a bit
      tags$br(),
      conditionalPanel(condition = "output.ifValidUploadedFile & input.computeLog2FC == 0",
                       div(textOutput('textLog2FC_2'), style = 'color:IndianRed;font-weight:bold;font-size:110%;text-align:center')),
      conditionalPanel(condition = "input.computeLog2FC",
                       div(style = "height: 100vh; width: 100%;",
                           plotlyOutput('plotNetwork') %>%
                             withSpinner(color="#56070C"),
                           div(style = "width: 100%; margin-top: 405px; display: flex; justify-content: center;",
                               downloadButton("downloadNetworkPlot", "Download network Plot")))
      )
    ),
    tabPanel(
      'Log',
      conditionalPanel(condition = "output.ifValidUploadedFile",
                       fluidRow(
                         column(width = 7, HTML('<br>
                                                <h3 style="margin-left: 2%;">Processing History</h3>
                                                <h5 style="margin-left: 2%;">This section logs processing
                                                parameters used in different trials (every time you hit Process button).</h5>
                                                <h5 style="margin-left: 2%;"><strong>Instructions:</strong> Click on a row to
                                                view the processed data at the specific point.</h5>
                                                <br>')),
                         column(width = 2, div(style = "width: 60%; margin-top: 30%;",
                                               actionButton('clearParamLog', 'Clear', width = '100%'))),
                         bsTooltip('clearParamLog', 'The parameter log table will be cleared.')
                       ),
                       fluidRow(
                         column(width = 7, DT::dataTableOutput("tblParamLog")),
                         column(width = 5,
                                conditionalPanel(condition = "typeof input.tblParamLog_rows_selected === 'undefined' || input.tblParamLog_rows_selected.length == 0",
                                                 HTML('<br>
                                                      <h4 style="color: IndianRed; text-align: left;">Please select a row.</h4>')
                                ),
                                conditionalPanel(condition = "typeof input.tblParamLog_rows_selected !== 'undefined' && input.tblParamLog_rows_selected.length > 0",
                                                 bsCollapse(
                                                   multiple = T,
                                                   bsCollapsePanel('Features Removed', style = 'info',
                                                                   textOutput('textRmFeatsLog') %>%
                                                                     withSpinner(color="#56070C")),
                                                   bsCollapsePanel('Data distribution', style = 'info',
                                                                   plotlyOutput('plotDatDistLog') %>%
                                                                     withSpinner(color="#56070C")),
                                                   bsCollapsePanel('Data completeness', style = 'info',
                                                                   plotlyOutput('plotDatCompleteLog') %>%
                                                                     withSpinner(color="#56070C")),
                                                   bsCollapsePanel('Quantification status', style = 'info',
                                                                   plotlyOutput('plotQuanStatusLog') %>%
                                                                     withSpinner(color="#56070C"))
                                                 )
                                )
                         )
                       )
      )
    )
  )
)

server <- function(input, output, session) {
  # Create reactive objects for storing up-to-date data
  reactMetabObj <- reactiveValues(metabObj = NULL, oriMetabObj = NULL, tmpMetabObj = NULL)
  reactLog2FCTbl <- reactiveVal()
  reactVulcanoHighlight <- reactiveVal()
  # Create reactive parameter list for detecting any parameter change by users,
  # which is to make data always follow built processing procedure
  reactParamList <- reactiveValues(smpFiltering = c(), featFiltering = c(), featCompleteCutoff = 0,
                                   featValidCutoff = 0, featValidStatus = c(), imputation = F,
                                   normalization = 'None')
  # Create reactive container to store processing parameters and processed data
  # of different trials for logging history and comparing results
  reactAnalysisLog <- reactiveValues(
    paramLogTbl = data.frame(idx = integer(),
                             smpFiltering = character(),
                             featClassFiltering = character(),
                             featFilteringNum = numeric(),
                             featCompleteCutoff = numeric(),
                             featValidCutoff = numeric(),
                             featValidStatus = character(),
                             imputation = character(),
                             normalization = character(),
                             stringsAsFactors = FALSE),
    rmFeatList = list(),
    # metabObjList = list(),
    # smpMetadatTblList = list(),
    metabAggreTblList = list()
  )
  
  # Show fileInput only when example data is not used
  output$updateFileInput <- renderUI({
    if (!input$exampleFile) {
      fileInput('uploadedFile', NULL, multiple = F, accept = '.xlsx',
                placeholder = 'No .xlsx file selected')
    }
  })
  # Show hint for uploading data when example data is being used
  output$textFileInput <- renderText({
    if (input$exampleFile) {
      'Uncheck example data box to upload your own data.'
    }
  })
  # Bring UI back to origin when example data box is unchecked
  observe({
    req(!input$exampleFile)
    reactMetabObj$metabObj <- NULL
  })
  # Initialize MetAlyzer SE object with example data
  observeEvent(input$exampleFile, {
    req(input$exampleFile)
    metabObj <- MetAlyzer_dataset(file_path = example_extraction_data(), silent = T)
    reactMetabObj$metabObj <- metabObj
    # Make copy of original data for filtering reset
    reactMetabObj$oriMetabObj <- metabObj
    # Make copy of original data for sanity check of filtering
    reactMetabObj$tmpMetabObj <- metabObj
    # Monitor whether example data is used to show further operations on client side
    output$ifValidUploadedFile <- reactive({
      !is.null(reactMetabObj$metabObj)
    })
    outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
    
    # Set parameters back to default
    updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
    updateSliderInput(session, 'featValidCutoffFiltering', value = 50)
    updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
    updateMaterialSwitch(session, 'imputation', value = T)
    updateSelectInput(session, 'normalization', selected = 'Median log₂ normalization')
    
    # Empty parameter log
    reactParamList$smpFiltering <- c()
    reactParamList$featFiltering <- c()
    reactParamList$featCompleteCutoff <- 0
    reactParamList$featValidCutoff <- 0
    reactParamList$featValidStatus <- c()
    reactParamList$imputation <- F
    reactParamList$normalization <- 'None'
  })
  # Initialize MetAlyzer SE object with uploaded data
  observeEvent(input$uploadedFile, {
    validUploadedFile <- try(
      metabObj <- MetAlyzer_dataset(file_path = input$uploadedFile$datapath,
                                    sheet = 1, silent = T),
      silent = T)
    if (!is(validUploadedFile, 'try-error')) {
      reactMetabObj$metabObj <- metabObj
      # Make copy of original data for filtering reset
      reactMetabObj$oriMetabObj <- metabObj
      # Make copy of original data for sanity check of filtering
      reactMetabObj$tmpMetabObj <- metabObj
      # Monitor whether uploaded file is valid to show further operations on client side
      output$ifValidUploadedFile <- reactive({
        !is.null(reactMetabObj$metabObj)
      })
      outputOptions(output, 'ifValidUploadedFile', suspendWhenHidden = F)
      
      # Set parameters back to default
      updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
      updateSliderInput(session, 'featValidCutoffFiltering', value = 50)
      updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
      updateMaterialSwitch(session, 'imputation', value = T)
      updateSelectInput(session, 'normalization', selected = 'Median log₂ normalization')
      
      # Empty parameter log
      reactParamList$smpFiltering <- c()
      reactParamList$featFiltering <- c()
      reactParamList$featCompleteCutoff <- 0
      reactParamList$featValidCutoff <- 0
      reactParamList$featValidStatus <- c()
      reactParamList$imputation <- F
      reactParamList$normalization <- 'None'
    } else {
      showModal(modalDialog(
        title = 'Uploaded file reading failed...',
        'Is the uploaded .xlsx file generated from MetIDQ™ software?',
        easyClose = T,
        footer = NULL
      ))
      reactMetabObj$metabObj <- NULL
    }
  })
  
  
  # Retrieve abundance data and sample metadata and compute feature completeness
  # level and quantification status validity level for showing data overviews
  datOverviewPack <- reactive({
    req(reactMetabObj$metabObj)
    metabAggreTbl <- metadata(reactMetabObj$metabObj)$aggregated_data %>%
      #### Do not know why Metabolite column is grouped by MetAlyzer
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
  
  
  # Prepare choices for feature filtering and log2(FC) vulcano highlighting
  featChoices <- reactive({
    req(reactMetabObj$metabObj)
    list(Class = as.list(unique(rowData(reactMetabObj$tmpMetabObj)$metabolic_classes)),
         Metabolite = as.list(rownames(reactMetabObj$tmpMetabObj)))
  })
  # Prepare needed information for sample filtering (i.e., choices, dictionary for
  # choice group searches, and list of choice groups containing duplicated choices)
  # and log2(FC) calculation
  smpChoicePack <- reactive({
    req(reactMetabObj$metabObj)
    smpChoiceList <- as.list(colData(reactMetabObj$tmpMetabObj))
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
    if (nrow(reactMetabObj$metabObj) == nrow(reactMetabObj$oriMetabObj)) {
      if ('Metabolism Indicators' %in% unlist(featChoices())) {
        updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
                          selected = 'Metabolism Indicators')
      } else {
        updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
      }
    }
  })
  # Update choices for sample filtering
  observe({
    req(smpChoicePack()$smpChoiceList)
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    if (ncol(reactMetabObj$metabObj) == ncol(reactMetabObj$oriMetabObj)) {
      updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
    }
  })
  
  # Create reactive values to monitor if data imputation and normalization are done
  # to avoid data normalized more than one time
  doneImputation <- reactiveVal(0)
  doneNormalization <- reactiveVal(0)
  # Create reactive values to monitor if any parameter is changed
  ifParamChange <- reactiveVal(0)
  
  # Rerun data processing if any parameter is changed, so that processing can follow
  # order from sample filtering, feature filtering, imputation, to normalization
  # Note that every first trial must be logged because output of selectInput will
  # never be c(), instead, NULL when nothing is selected
  observeEvent(input$updateProcessing, {
    if (any(!identical(sort(reactParamList$smpFiltering), sort(input$smpChoicesFiltering)),
            !identical(sort(reactParamList$featFiltering), sort(input$featChoicesFiltering)),
            !identical(reactParamList$featCompleteCutoff, input$featCompleteCutoffFiltering),
            !identical(reactParamList$featValidCutoff, input$featValidCutoffFiltering),
            !identical(reactParamList$featValidStatus, input$featValidStatusFiltering),
            !identical(reactParamList$imputation, input$imputation),
            !identical(reactParamList$normalization, input$normalization))) {
      # Revert processed temporary MetAlyzer object and parameter log back to origins
      reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
      doneImputation(0)
      doneNormalization(0)
      
      reactParamList$smpFiltering <- c()
      reactParamList$featFiltering <- c()
      reactParamList$featCompleteCutoff <- 0
      reactParamList$featValidCutoff <- 0
      reactParamList$featValidStatus <- c()
      reactParamList$imputation <- F
      reactParamList$normalization <- 'None'
      
      # Turn reactive ifParamChange to 1, so used parameters and processed data
      # will be logged
      ifParamChange(1)
    }
  })
  # Do sample filtering (place this chunk before 'Do feature filtering', so that
  # feature filtering is executed based on sample filtered data)
  observeEvent(input$updateProcessing, {
    req(input$smpChoicesFiltering)
    for (selectedChoice in input$smpChoicesFiltering) {
      selectedChoiceGp <- smpChoicePack()$smpChoices2Gps[selectedChoice]
      # Prepare up-to-date choices of selected choice group to avoid error that
      # metadata of same sample is selected multiple times at once
      choices <- unique(colData(reactMetabObj$tmpMetabObj)[[selectedChoiceGp]])
      # Convert NA into character so that it can be searched
      if (any(is.na(choices))) {
        choices <- choices[-which(is.na(choices))]
        choices <- c(choices, 'NA')
      }
      if (!selectedChoiceGp %in% smpChoicePack()$dupSmpChoiceGps) {
        if (selectedChoice %in% choices) {
          if (selectedChoice == 'NA') {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetaData(reactMetabObj$tmpMetabObj,
                                                                   !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetaData(reactMetabObj$tmpMetabObj,
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
            reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetaData(reactMetabObj$tmpMetabObj,
                                                                   !is.na(.data[[selectedChoiceGp]]))
          } else {
            reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetaData(reactMetabObj$tmpMetabObj,
                                                                   !.data[[selectedChoiceGp]] %in% modSelectedChoice)
          }
        }
      }
    }
    # Avoid app crash when no sample is left
    if (ncol(reactMetabObj$tmpMetabObj) == 0) {
      showModal(modalDialog(
        title = 'No sample is left after sample filtering...',
        'Please respecify parameters for sample filtering.',
        easyClose = T,
        footer = NULL
      ))
    }
  })
  # Do feature filtering
  #### Add advanced feature filtering: Modified 80% rule
  observeEvent(input$updateProcessing, {
    req(datOverviewPack())
    # Collect features to remove based on user's selection
    rmSelectedFeats <- input$featChoicesFiltering
    # Collect features to remove based on missingness
    featCompleteLvTbl <- datOverviewPack()$featCompleteLvTbl
    featCompleteLevels <- featCompleteLvTbl$CompleteRatio
    featCompleteCutoff <- input$featCompleteCutoffFiltering / 100
    rmMissFeats <- featCompleteLvTbl$Metabolite[featCompleteLevels < featCompleteCutoff]
    # Collect features to remove based on validity
    # featStatusCountTbl <- datOverviewPack()$featStatusCountTbl
    featValidCutoff <- input$featValidCutoffFiltering / 100
    featValidStatus <- input$featValidStatusFiltering
    # featValidCounts <- rep(0, nrow(featStatusCountTbl))
    # for (status in featValidStatus) {
    #   featValidCounts <- featValidCounts + featStatusCountTbl[[paste0(status, 'Count')]]
    # }
    # featValidLevels <- featValidCounts / featStatusCountTbl[['TotalCount']]
    # rmInvalidFeats <- featStatusCountTbl$Metabolite[featValidLevels < featValidCutoff]
    # # Summarize features to remove
    # rmFeats <- unique(c(rmSelectedFeats, rmMissFeats)) # rmInvalidFeats
    # Skip feature filtering if no sample was left in sample filtering part since
    # feature filtering is based on sample filtered data
    if (ncol(reactMetabObj$tmpMetabObj) != 0) {
      #### Filtering results depend on whether rmSelectedFeats is Null when specifying
      #### rmFeats to 'drop_metabolites'????
      reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetabolites(reactMetabObj$tmpMetabObj,
                                                                drop_metabolites = rmSelectedFeats,
                                                                drop_NA_concentration = FALSE)
      reactMetabObj$tmpMetabObj <- MetAlyzer::filterMetabolites(reactMetabObj$tmpMetabObj,
                                                                drop_metabolites = rmMissFeats,
                                                                drop_NA_concentration = FALSE,
                                                                min_percent_valid = featValidCutoff,
                                                                valid_status = featValidStatus,
                                                                per_group = NULL)
      # Avoid app crash when no feature is left, e.g., min_percent_valid > 0 and valid_status == c()
      if (nrow(reactMetabObj$tmpMetabObj) == 0) {
        showModal(modalDialog(
          title = 'No metabolite is left after feature filtering...',
          'Please respecify parameters for metabolite filtering.',
          easyClose = T,
          footer = NULL
        ))
      }
    }
  })
  # Do data imputation and normalization
  observeEvent(input$updateProcessing, {
    req(reactMetabObj$tmpMetabObj)
    # Skip imputation and normalization if no sample or feature was left after filtering
    if (all(ncol(reactMetabObj$tmpMetabObj) != 0, nrow(reactMetabObj$tmpMetabObj) != 0)) {
      if (all(input$imputation, doneImputation() == 0)) {
        reactMetabObj$tmpMetabObj <- data_imputation(reactMetabObj$tmpMetabObj)
        doneImputation(1)
      }
      if (doneNormalization() == 0) {
        if (input$normalization == 'Total ion count (TIC) log₂ normalization') {
          reactMetabObj$tmpMetabObj <- data_normalization(reactMetabObj$tmpMetabObj, norm_method = 'TIC')
          doneNormalization(1)
        } else if (input$normalization == 'Median log₂ normalization') {
          reactMetabObj$tmpMetabObj <- data_normalization(reactMetabObj$tmpMetabObj, norm_method = 'median')
          doneNormalization(1)
        }
      }
      # Update main MetAlyzer object for further analysis
      reactMetabObj$metabObj <- reactMetabObj$tmpMetabObj
      
      
      # Log parameters and processed data
      if (ifParamChange() == 1) {
        # Log parameters
        reactParamList$smpFiltering <- input$smpChoicesFiltering
        reactParamList$featFiltering <- input$featChoicesFiltering
        reactParamList$featCompleteCutoff <- input$featCompleteCutoffFiltering
        reactParamList$featValidCutoff <- input$featValidCutoffFiltering
        reactParamList$featValidStatus <- input$featValidStatusFiltering
        reactParamList$imputation <- input$imputation
        reactParamList$normalization <- input$normalization
        
        
        # Prepare removed feature list
        selectedFeats <- input$featChoicesFiltering
        featClasses <- unique(rowData(reactMetabObj$oriMetabObj)$metabolic_classes)
        rmFeatClasses <- selectedFeats[selectedFeats %in% featClasses]
        # Prepare removed feature class list
        oriFeatSpace <- rownames(reactMetabObj$oriMetabObj)
        filtFeatSpace <- rownames(reactMetabObj$metabObj)
        rmFeats <- oriFeatSpace[!oriFeatSpace %in% filtFeatSpace]
        
        # Prevent empty cells for better readability
        selectedSmps <- input$smpChoicesFiltering
        if (length(selectedSmps) == 0) {
          selectedSmps <- 'None'
        }
        
        if (length(rmFeatClasses) == 0) {
          rmFeatClasses <- 'None'
        }
        
        selectedValidStatus <- input$featValidStatusFiltering
        if (length(selectedValidStatus) == 0) {
          selectedValidStatus <- 'None'
        }
        
        newParamLog <- data.frame(
          idx = nrow(reactAnalysisLog$paramLogTbl) + 1,
          #### What is effect of I(list())?
          smpFiltering = I(list(paste(selectedSmps, collapse = ' | '))),
          featClassFiltering = paste(rmFeatClasses, collapse = ' | '),
          featFilteringNum = length(rmFeats),
          featCompleteCutoff = reactParamList$featCompleteCutoff,
          featValidCutoff = reactParamList$featValidCutoff,
          featValidStatus = I(list(paste(selectedValidStatus, collapse = ' | '))),
          imputation = as.character(reactParamList$imputation),
          normalization = reactParamList$normalization,
          stringsAsFactors = FALSE
        )
        reactAnalysisLog$paramLogTbl <- rbind(reactAnalysisLog$paramLogTbl, newParamLog)
        
        # Log processed data
        logIdx <- nrow(reactAnalysisLog$paramLogTbl)
        reactAnalysisLog$rmFeatList[[logIdx]] <- rmFeats
        # reactAnalysisLog$metabObjList[[logIdx]] <- reactMetabObj$metabObj
        # reactAnalysisLog$smpMetadatTblList[[logIdx]] <- datOverviewPack()$smpMetadatTbl
        reactAnalysisLog$metabAggreTblList[[logIdx]] <- datOverviewPack()$metabAggreTbl
        
        # Return reactive ifParamChange to 0, so this chunk will run only if any
        # parameter is changed
        ifParamChange(0)
      } 
    } else {
      # Revert temporary MetAlyzer object to origin for redoing filtering
      reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
    }
  })
  # Revert processed data and specified parameters back to origins
  observeEvent(input$revertProcessing, {
    req(reactMetabObj$metabObj)
    reactMetabObj$metabObj <- reactMetabObj$oriMetabObj
    reactMetabObj$tmpMetabObj <- reactMetabObj$oriMetabObj
    doneImputation(0)
    doneNormalization(0)
    ifParamChange(0)
    
    # Set parameters for feature filtering back to default
    if ('Metabolism Indicators' %in% unlist(featChoices())) {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices(),
                        selected = 'Metabolism Indicators')
    } else {
      updateSelectInput(session, 'featChoicesFiltering', choices = featChoices())
    }
    updateSliderInput(session, 'featCompleteCutoffFiltering', value = 80)
    updateSliderInput(session, 'featValidCutoffFiltering', value = 50)
    updateCheckboxGroupInput(session, 'featValidStatusFiltering', selected = c('Valid', 'LOQ'))
    # Set parameters for sample filtering back to default
    # Make choices lists so that sole choice in certain choice group can be shown
    smpChoiceList <- lapply(smpChoicePack()$smpChoiceList, function(choices) {
      as.list(choices)
    })
    updateSelectInput(session, 'smpChoicesFiltering', choices = smpChoiceList)
    # Set parameters for imputation and normalization back to default
    updateMaterialSwitch(session, 'imputation', value = T)
    updateSelectInput(session, 'normalization', selected = 'Median log₂ normalization')
    
    # Empty parameter log
    reactParamList$smpFiltering <- c()
    reactParamList$featFiltering <- c()
    reactParamList$featCompleteCutoff <- 0
    reactParamList$featValidCutoff <- 0
    reactParamList$featValidStatus <- c()
    reactParamList$imputation <- F
    reactParamList$normalization <- 'None'
  })
  # Clear parameter log table
  observeEvent(input$clearParamLog, {
    reactAnalysisLog$paramLogTbl <- data.frame(idx = integer(),
                                               smpFiltering = character(),
                                               featClassFiltering = character(),
                                               featFilteringNum = numeric(),
                                               featCompleteCutoff = numeric(),
                                               featValidCutoff = numeric(),
                                               featValidStatus = character(),
                                               imputation = character(),
                                               normalization = character(),
                                               stringsAsFactors = FALSE)
    reactAnalysisLog$rmFeatList <- list()
    # reactAnalysisLog$metabObjList = list()
    # reactAnalysisLog$smpMetadatTblList = list()
    reactAnalysisLog$metabAggreTblList = list()
    # Turn reactive ifParamChange to 1, so current unchanged parameters can be logged
    # if Process button is hit
    ifParamChange(1)
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
      metabObj <- reactMetabObj$metabObj
      # Do log2 transformation on data for 'calc_log2FC' if normalization was not performed
      if (doneNormalization() == 0) {
        oriConc <- metadata(metabObj)$aggregated_data$Concentration
        metadata(metabObj)$aggregated_data$Concentration <- glog2(oriConc)
      }
      # Extract samples of interest
      selectedChoiceGp <- input$smpChoiceGpsLog2FC
      selectedChoices <- c(input$smpChoicesLog2FC_1, input$smpChoicesLog2FC_2)
      if (!'NA' %in% selectedChoices) {
        metabObj <- MetAlyzer::filterMetaData(metabObj, .data[[selectedChoiceGp]] %in% selectedChoices)
      } else {
        metabObj <- MetAlyzer::filterMetaData(metabObj, is.na(.data[[selectedChoiceGp]]) |
                                                .data[[selectedChoiceGp]] %in% selectedChoices)
      }
      metabObj <- calc_log2FC(metalyzer_se = metabObj,
                              categorical = selectedChoiceGp)
      reactLog2FCTbl(log2FC(metabObj))
      
      # Update the slider input, for custom inputs
      updateSliderInput(session, "plotVolcanoLog2FCCutoff",
                        max = floor(max(na.omit(reactLog2FCTbl()$log2FC))))
    } else {
      showModal(modalDialog(
        title = 'Log₂(FC) computation failed...',
        'Please select two different sample groups.',
        easyClose = T,
        footer = NULL
      ))
    }
  })
  
  
  # Update choices for sample metadata to include in exported abundance csv file
  observe({
    req(reactMetabObj$metabObj)
    updateSelectInput(session, 'metaChoicesRawAbunExport', choices = colnames(colData(reactMetabObj$metabObj)))
  })
  # Update choices for feature identifiers to include in exported annotation csv file
  featIdTbl <- reactive({
    readxl::read_excel('data/mappingTable.xlsx')
  })
  output$loadFeatIdChoicesExport <- renderUI({
    annoChoices <- colnames(featIdTbl())
    # Remove columns that will be included in exported file
    annoChoices <- annoChoices[-which(annoChoices %in% c('BiocratesName', 'Match'))]
    selectInput('featIdChoicesExport', 'Select metabolite annotation(s) to include:',
                choices = annoChoices, multiple = T)
  })
  
  
  # Set plot theme
  # th <- theme_bw(base_size = 15) +
  #   theme(axis.title = element_text(face = 'bold'),
  #         axis.text = element_text(face = 'bold'),
  #         axis.ticks = element_line(linewidth = 0.8),
  #         legend.text = element_text(size = 15))
  
  # Display Log Elements
  # Log of parameters
  output$tblParamLog <- DT::renderDataTable({
    paramLogTbl <- reactAnalysisLog$paramLogTbl
    # Reverse parameter log table to display latest log on top
    paramLogTbl <- paramLogTbl[order(-paramLogTbl$idx),] %>%
      dplyr::mutate(imputation = dplyr::case_when(imputation %in% 'TRUE' ~ 'True',
                                                  imputation %in% 'FALSE' ~ 'False'),
                    featCompleteCutoff = paste0(featCompleteCutoff, '%'),
                    featValidCutoff = paste0(featValidCutoff, '%'))
    colnames(paramLogTbl) <- c('Process', 'Samples Removed', 'Feature Classes Removed',
                               'Count of Features Removed', 'Completeness Cutoff',
                               'Validity Cutoff', 'Validity Statuses', 'Imputation', 'Normalization')
    DT::datatable(paramLogTbl, rownames = F, options = list(dom = 't', ordering = F),
                  selection = list(mode = 'single', target = 'row'), style = 'bootstrap')
  })
  
  # Prepare row index clicked
  clickedRowIdx <- reactive({
    req(input$tblParamLog_rows_selected)
    clickedRow <- input$tblParamLog_rows_selected
    # Because of reversed parameter log table
    sort(reactAnalysisLog$paramLogTbl$idx, decreasing = T)[clickedRow]
  })
  
  # Log of all removed features
  output$textRmFeatsLog <- renderText({
    req(clickedRowIdx())
    rmFeats <- reactAnalysisLog$rmFeatList[[clickedRowIdx()]]
    # Place it here so that empty removed feature list got length of 0 in parameter log table
    if (length(rmFeats) == 0) {
      rmFeats <- 'None'
    }
    paste(sort(rmFeats), collapse = ' | ')
  })
  
  # Log of processed data distribution
  output$plotDatDistLog <- renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    g <- ggplot(actMetabAggreTbl, aes(x=ID, y=Concentration)) +
      geom_boxplot() +
      labs(x = 'Sample', y = 'Metabolite abundance') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Put abundance on log10 scale for better visualization if data is not normalized 
    if (reactAnalysisLog$paramLogTbl$normalization[clickedRowIdx()] == 'None') {
      g <- g + scale_y_log10()
    }
    ggplotly(g)
  })
  
  # Log of processed data quantification status
  output$plotQuanStatusLog <- renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    smpStatusCountTbl <- actMetabAggreTbl %>%
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
  
  # Log of processed data completeness
  output$plotDatCompleteLog <- renderPlotly({
    req(clickedRowIdx())
    actMetabAggreTbl <- reactAnalysisLog$metabAggreTblList[[clickedRowIdx()]]
    smpCompleteCountTbl <- actMetabAggreTbl %>%
      dplyr::mutate(Completeness = dplyr::case_when(!Concentration %in% c(0, NA) ~ 'Observed',
                                                    Concentration %in% c(0, NA) ~ 'Missing')) %>%
      dplyr::group_by(ID, Completeness) %>%
      dplyr::summarise(Count = dplyr::n()) %>%
      dplyr::ungroup()
    ggplotly(
      ggplot(smpCompleteCountTbl, aes(x=ID, y=Count, fill=Completeness)) +
        geom_col(position = 'stack') +
        scale_fill_manual(values = c(Missing = 'grey', Observed = 'black')) +
        labs(x = 'Sample') +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )
  })
  
  
  # Show data overviews
  # Data distribution
  output$updateGpColsDatDist <- renderUI({
    req(smpChoicePack()$smpChoiceList)
    smpChoiceList <- smpChoicePack()$smpChoiceList
    # Remove choice groups whose level sizes are one or equal to sample size
    choiceGpSizes <- sapply(smpChoiceList, length)
    rmChoiceGps <- which(choiceGpSizes == 1 | choiceGpSizes == ncol(reactMetabObj$metabObj))
    if (length(rmChoiceGps) != 0) {
      smpChoiceGps <- names(smpChoiceList)[-rmChoiceGps]
    } else {
      smpChoiceGps <- names(smpChoiceList)
    }
    selectInput('gpColsDatDist', 'Color by:',
                choices = c('None', smpChoiceGps),
                selected = 'None', multiple = F)
  })
  output$plotDatDist <- renderPlotly({
    req(datOverviewPack()$metabAggreTbl, input$gpColsDatDist)
    metabAggreTbl <- datOverviewPack()$metabAggreTbl
    if (input$gpColsDatDist == 'None') {
      g <- ggplot(metabAggreTbl, aes(x=ID, y=Concentration)) +
        geom_boxplot()
    } else {
      g <- ggplot(metabAggreTbl, aes(x=ID, y=Concentration, fill=.data[[input$gpColsDatDist]])) +
        geom_boxplot(alpha = 1) +
        scale_fill_brewer(palette = 'Set1')
    }
    g <- g +
      labs(x = 'Sample', y = 'Metabolite abundance') +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Put abundance on log10 scale for better visualization if data is not normalized 
    if (doneNormalization() == 0) {
      g <- g + scale_y_log10()
    }
    ggplotly(g)
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
        scale_fill_manual(values = c(Missing = 'grey', Observed = 'black')) +
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
    paste(sum(featStatusValid$ValidRatio < 0.5), 'out of', nrow(reactMetabObj$metabObj),
          'metabolites have less than 50% measurements with valid status (Valid, LOQ),',
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
  
  # Export abundance matrix
  output$downloadRawAbunExport <- downloadHandler(
    filename = paste("abundance_matrix_", format(Sys.Date(), "%Y-%m-%d"), ".csv", sep = ""),
    content = function(file) {
      #### Add option of exporting also processed data?
      #### What inputs MetaboAnalyst needs?
      MetAlyzer::exportAbunValues(reactMetabObj$metabObj,
                                  input$metaChoicesRawAbunExport,
                                  file_path = file)
    }
  )
  # Export metabolite identifiers
  output$downloadFeatIdsExport <- downloadHandler(
    filename = paste("metabolite_identifiers_", format(Sys.Date(), "%Y-%m-%d"), ".csv", sep = ""),
    content = function(file) {
      exportfeatIdTbl <- featIdTbl()[featIdTbl()$BiocratesName %in% rownames(reactMetabObj$metabObj),]
      exportfeatIdTbl <- dplyr::select(exportfeatIdTbl, BiocratesName, Match, input$featIdChoicesExport)
      write.csv(exportfeatIdTbl, file)
    }
  )
  
  
  # Visualize log2(FC)
  # Give sign before log2(FC) calculation
  output$textLog2FC <- renderText({
    'Please calculate Log₂(FC) first.'
  })
  output$textLog2FC_2 <- renderText({
    'Please calculate Log₂(FC) first.'
  })
  
  # Update choices for metabolite highlighting
  observe({
    req(featChoices())
    updateSelectInput(session, 'metabChoicesVulcano', choices = featChoices())
  })
  # Create new column for highlighting metabolites in vulcano and scatter plots
  observe({
    req(reactLog2FCTbl())
    selectedChoices <- input$metabChoicesVulcano
    highlightTbl <- reactLog2FCTbl()
    if (!is.null(selectedChoices)) {
      highlightTbl$highlight <- FALSE
      highlightTbl$highlight[highlightTbl$Metabolite %in% selectedChoices] <- TRUE
      highlightTbl$highlight[highlightTbl$Class %in% selectedChoices] <- TRUE
    } else if (is.null(selectedChoices) | input$highlightVulcano) {
      highlightTbl$highlight <- FALSE
    }
    reactVulcanoHighlight(highlightTbl)
  })
  
  # Vulcano plot
  output$plotVolcano <- renderPlotly({
    req(reactLog2FCTbl())
    if (!input$highlightVulcano) {
      plotly_vulcano(reactLog2FCTbl(),
                     cutoff_x = input$plotVolcanoLog2FCCutoff,
                     cutoff_y = as.numeric(input$plotVolcanoPValCutoff))
    } else {
      plotly_vulcano(reactVulcanoHighlight(),
                     cutoff_x = input$plotVolcanoLog2FCCutoff,
                     cutoff_y = as.numeric(input$plotVolcanoPValCutoff))
    }
  })
  
  # Scatter plot
  output$plotScatter <- renderPlotly({
    req(reactLog2FCTbl())
    # if (!input$highlightVulcano) {
    #   plot <- plotly_scatter(reactLog2FCTbl())
    # } else {
    #   plot <- plotly_scatter(reactVulcanoHighlight())
    # }
    plot <- plotly_scatter(reactLog2FCTbl())
    hide_legend(plot$Plot)
  })
  output$plotScatterLegend <- renderImage({
    req(reactLog2FCTbl()) 
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
         alt = "",
         deleteFile = TRUE)
  }, deleteFile = TRUE)
  
  # Network plot
  output$plotNetwork <- renderPlotly({
    req(reactLog2FCTbl())
    plotly_network(reactLog2FCTbl())
  })
  
  # Download the log2(FC) visuals as html
  output$downloadVulcanoPlot <- downloadHandler(
    filename = function() {
      "vulcano_plot.html"
    },
    content = function(file) {
      # Define a variable to store the final plot
      if (input$highlightVulcano) {
        req(reactVulcanoHighlight())
        final_plot <- plotly_vulcano(reactVulcanoHighlight(), 
                                     cutoff_x = input$plotVolcanoLog2FCCutoff,
                                     cutoff_y = as.numeric(input$plotVolcanoPValCutoff))
      } else {
        final_plot <- plotly_vulcano(reactLog2FCTbl(), 
                                     cutoff_x = input$plotVolcanoLog2FCCutoff,
                                     cutoff_y = as.numeric(input$plotVolcanoPValCutoff))
      }
      
      # Save the Plotly plot as an HTML file
      htmlwidgets::saveWidget(
        widget = final_plot,
        file = file,
        selfcontained = TRUE
      )
    }
  )
  output$downloadScatterPlot <- downloadHandler(
    filename = function() {
      "scatter_plot.html"
    },
    content = function(file) {
      # Define a variable to store the final plot
      # if (input$highlightVulcano) {
      #   req(reactVulcanoHighlight())
      #   final_plot <- plotly_scatter(reactVulcanoHighlight())$Plot
      # } else {
      #   final_plot <- plotly_scatter(reactLog2FCTbl())$Plot
      # }
      final_plot <- plotly_scatter(reactLog2FCTbl())$Plot
      
      # Save the Plotly plot as an HTML file
      htmlwidgets::saveWidget(
        widget = final_plot,
        file = file,
        selfcontained = TRUE
      )
    }
  )
  output$downloadNetworkPlot <- downloadHandler(
    filename = function() {
      "network_plot.html"
    },
    content = function(file) {
      # Save the Plotly plot as an HTML file
      htmlwidgets::saveWidget(
        widget = plotly_network(reactLog2FCTbl()),
        file = file,
        selfcontained = TRUE
      )
    }
  )
}

shinyApp(ui = ui, server = server)
